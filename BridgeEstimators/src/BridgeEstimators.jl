module BridgeEstimators

export BirthDeathProcess, OrnsteinUhlenbeckProcess
export sample, subsample, likelihood

using SimpleBirthDeathProcess
using Distributions, Random


### Observations
abstract type Observations{T} end
struct RegularObservations{T} <: Observations{T}
    values::Array{T}
    δt::Real
end
struct IrregularObservations{T} <: Observations{T}
    values::Array{T}
    times::Vector{Real}
end



subsample(X::RegularObservations, i::Int) = RegularObservations(X.values[:, i], X.δt)
subsample(X::IrregularObservations, i::Int) = RegularObservations(X.values[:, i], X.times)

bitVector = Union{Vector{Bool}, BitVector}
subsample(X::RegularObservations, I::bitVector) = RegularObservations(X.values[:, I], X.δt)
subsample(X::IrregularObservations, I::bitVector) = RegularObservations(X.values[:, I], X.times)


### Processes
abstract type Process{T} end
struct BirthDeathProcess <: Process{Int}
    parameters::Vector # [λ, μ]
    initial::Int
end
struct OrnsteinUhlenbeckProcess <: Process{Real}
    parameters::Vector # [θ, σ]
    initial::Real
end

sample(P::BirthDeathProcess, δt::Real, n::Int) = RegularObservations{Int}(rand_discrete(n, P.initial, n, δt, P.parameters).state, δt)
sample(npaths::Int, P::BirthDeathProcess, δt::Real, nobs::Int) = RegularObservations{Int}(rand_discrete(npaths, P.initial, nobs, δt, P.parameters).state, δt)

function likelihood(P::BirthDeathProcess, X::RegularObservations; log = true, bridge = false) 
    l = loglik(P.parameters, observation_discrete_time_even(X.δt, X.values)
)
    if bridge 
        l -= trans_prob(X.values[1], X.values[end], X.δt * (size(X.values, 1) - 1), P.parameters)
    end

    return log ? l : exp(l)
end


function kernel(P::OrnsteinUhlenbeckProcess, x, δt) 
    β, σ = P.parameters
    return Normal(x*exp(-β*δt), sqrt((0.5σ^2/β) *(1-exp(-2*β*δt))))
end

function sample(P::OrnsteinUhlenbeckProcess, δt::Real, n::Int)
    values = Vector{Float64}(undef, n)
    values[1] = P.initial
    for i in 1:n-1
        values[i+1] = rand(kernel(P, values[i], δt))
    end
    return RegularObservations(values, δt)
end

function sample(npaths::Int, P::OrnsteinUhlenbeckProcess, δt::Real, n::Int, seed = 1)
    Random.seed!(seed)
    values = Array{Float64}(undef, n, npaths)
    values[1, :] .= fill(P.initial, npaths)
    for i in 1:n-1, j in 1:npaths
        values[i+1, j] = rand(kernel(P, values[i, j], δt))
    end
    return RegularObservations(values, δt)
end

function likelihood(P::OrnsteinUhlenbeckProcess, X::RegularObservations; log = true, bridge = false) 
    n = size(X.values, 1)
    L = [logpdf(kernel(P, X.values[i], X.δt), X.values[i+1]) for i in 1:n-1]
    l = sum(L)
    if bridge 
        l -= logpdf(kernel(P, X.values[1], X.δt * (n-1)), X.values[end])
    end

    return log ? l : exp(l)
end



end # module
