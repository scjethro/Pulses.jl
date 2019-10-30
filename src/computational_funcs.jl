module computational_funcs
# some computational functions for calculating stuff

using DifferentialEquations:ODEProblem, solve, Tsit5
using LinearAlgebra:Matrix, I, mul!
using StaticArrays

export propagator, full_propagator, fast_propagator, full_fast_propagator, expectation_value, track_expectation_value

# in place von Neumann equation
function comp_U!(du, u, p, t)
    mul!(du, -1.0im * p.func(p, t), u)
end

# von Neumann equation for use with static arrays
function comp_U(u, p, t)
    -1.0im * p.func(p, t) * u
end

# computes the closed system propagator U
function propagator(p, td, u0::T)::T where T
    # p is a named tuple containing at least a field called func which computes the Hamiltonian
    tspan = (0.0, td)
    prob = ODEProblem(comp_U, u0, tspan, p)
    sol = solve(prob, abstol = 1e-6, reltol = 1e-6, save_everystep=false, alg = Tsit5())
    return sol.u[end]
end

# computes the propagator but returns the ODE solution object
function full_propagator(p, td, u0)
    tspan = (0.0, td)
    prob = ODEProblem(comp_U, u0, tspan, p)
    sol = solve(prob, abstol = 1e-6, reltol = 1e-6, save_everystep=true, alg = Tsit5())
    return sol
end

# this solves the ODE much faster *but* it doesn't specify the tolerances so you might see some deviation
function fast_propagator(p, td, u0::T)::T where T
    tspan = (0.0, td)
    prob = ODEProblem(comp_U, u0, tspan, p)
    sol = solve(prob, save_everystep=false, alg = Tsit5())
    return sol.u[end]
end

function full_fast_propagator(p, td, u0)
    tspan = (0.0, td)
    prob = ODEProblem(comp_U, u0, tspan, p)
    sol = solve(prob, save_everystep=true, alg = Tsit5())
    return sol
end

# computes the expectation value of an operator and a state psi
function expectation_value(psi, operator)::AbstractFloat
    real(psi' * operator * psi)
end

# track the expectation value of an operator as a function of time
function track_expectation_value(sol, psi, operator, res = 1000)
    z = @MVector zeros(res)
    for i = 1:res
        t = (i-1)/res * sol.t[end]
        @inbounds z[i] = expectation_value(sol(t) * psi, operator)
    end
    return z
end

end
