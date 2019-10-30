module Pulses

# more coming soon, once I rewrite them
include("computational_funcs.jl")
include("pulse_defs.jl")

export propagator, full_propagator, fast_propagator, full_fast_propagator,
expectation_value, track_expectation_value,
pulse, pi_pulse, pi_half_pulse, delay, optimal_pulse, PulseSequence, nmr_pulse, array_pulse

end
