module pulse_defs

using DelimitedFiles:readdlm
export pulse, pi_pulse, pi_half_pulse, delay, optimal_pulse, PulseSequence, nmr_pulse,
array_pulse

# everything is a subttype of Pulse
abstract type pulse end

# struct containing the data for a pi pulse
struct pi_pulse{T<:AbstractFloat} <: pulse
	amplitude::T
	duration::T
	phase::T

	function pi_pulse(amp::T, ph::T) where T<: AbstractFloat
		new{T}(amp, pi/amp, deg2rad(ph))
	end
end

struct pi_half_pulse{T<:AbstractFloat} <: pulse
	amplitude::T
	duration::T
	phase::T

	function pi_half_pulse(amp::T, ph::T) where T<:AbstractFloat
		new{T}(amp, pi/amp/2, deg2rad(ph))
	end
end

struct delay{T<:AbstractFloat} <: pulse
	amplitude::T
	duration::T
	phase::T

	function delay(dur::T) where T<:AbstractFloat
		new{T}(0.0, dur, 0.0)
	end
end

struct nmr_pulse{T<:AbstractFloat} <: pulse
	amplitude::T
	duration::T
	phase::T
	function nmr_pulse(amp::T, dur::T, ph::T) where T <: AbstractFloat
		new{T}(amp, dur, ph)
	end
end

# object that stores the data for a shaped pulse
struct optimal_pulse{S<:String, T<:AbstractFloat} <: pulse
	filename::S
	amplitude::Array{T,1}
	duration::T
	phase::Array{T,1}
	time::Array{T,1} # we need time to be linear here

	function optimal_pulse(filename::String, duration::AbstractFloat)
		f = open(filename, "r")
		m = readdlm(f, ' ')
		fake_time = collect(range(0.0, 1.0, length=length(m[:,2]))) * duration
		new{String, AbstractFloat}(filename, m[:,2] * 2.0 * pi, duration, m[:,3], fake_time)
	end

end

# a pulse with an analytical formula for phase + amplitude (or just an array for each)
struct array_pulse{T <: AbstractFloat} <: pulse
	amplitude::Array{T, 1}
	duration::T
	phase::Array{T, 1}
	time::Array{T, 1}

	function array_pulse(amp_arr::Array{T, 1}, duration::T, ph_arr::Array{T, 1}) where T <: AbstractFloat
		new{AbstractFloat}(amp_arr, duration, ph_arr, collect(range(0, duration, length = length(amp_arr))))
	end

	function array_pulse(amp_func, duration, ph_func, res::Integer)
		fake_time = range(0, 1, length = res) # fake time will be converted to the duration of the pulse
		new{AbstractFloat}(amp_func.(fake_time), duration, ph_func.(fake_time), fake_time .* duration)
	end
end

# stores a list of pulses but hopefully in a convenient way that speeds up simulations
struct PulseSequence{T<:Float64, P<:pulse}
	sequence::Array{P}
	rabi::T
	pulse_timings::Array{T, 1}
	cumulative_timings::Array{T, 1}
	total_duration::T

	function PulseSequence(seq::Array{pulse}, rab::Float64)
		pt = _get_pulse_timings(seq)
		new{Float64, pulse}(seq, rab, pt, cumsum(pt), sum(pt))
	end

	function _get_pulse_timings(seq)
		pt = []
		for (i, pulse) in enumerate(seq)
			append!(pt, pulse.duration)
		end
		return pt
	end
end

# needs to be reworked
function _amp_ph_t(c::PulseSequence, t::AbstractFloat)::Tuple{AbstractFloat,AbstractFloat}
	# it might be faster to just accurately type the return type of the slice?
	# I think the issue is that we aren't specifying the dimension of the view so we get typed with Any
	ind = searchsortedfirst(c.cumulative_timings, t)
	if ind > length(c.cumulative_timings)
		ind = length(c.cumulative_timings)
	end
	if ind >= 2
		# we take time from the beginning of the pulse
		t = t - c.cumulative_timings[ind - 1]
	end
	return c.sequence[ind](t)
end

function _amp_ph_t(c::Union{optimal_pulse, array_pulse}, t::AbstractFloat)::Tuple{AbstractFloat, AbstractFloat}
	(find_val_in_array(c.amplitude, t, c.duration), find_val_in_array(c.phase, t, c.duration))
end

function _amp_ph_t(c::pulse, t::AbstractFloat)::Tuple{AbstractFloat, AbstractFloat}
	(c.amplitude, c.phase)
end

# needs to be reworked
function find_val_in_array(array, t, t_stop)
	i = find_index_in_array(array, t, t_stop)
	return array[i]
end

function find_index_in_array(array, t, t_stop)
	i = Integer(floor(t/t_stop * length(array)))
    if i > length(array)
        return length(array)
    elseif i == 0
        return 1
    end
	return i
end

# lord forgive me for the things I do
(c::PulseSequence)(b::Float64) = _amp_ph_t(c, b)
(c::pi_pulse)(b::Float64) = _amp_ph_t(c, b)
(c::pi_half_pulse)(b::Float64) = _amp_ph_t(c, b)
(c::optimal_pulse)(b::Float64) = _amp_ph_t(c, b)
(c::array_pulse)(b::Float64) = _amp_ph_t(c, b)
(c::delay)(b::Float64) = _amp_ph_t(c, b)

end
