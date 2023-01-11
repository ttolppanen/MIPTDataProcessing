# using QuantumOperators
# using QuantumTimeEvolution
# using ITensors

function get_measurements()
    measurements = Dict(
        :std_n => (; exact = std_n_exact, mps = std_n_mps)
    )
    return measurements
end

function std_n_exact(; sp...)
    msrop = measurementoperators(nop(sp[:d]), sp[:L])
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob)
    return meffect!
end
function std_n_mps(ITensors_apply_kwargs; sp...)
    msrop = measurementoperators(nop(sp[:d]), siteinds(sp[:state0]))
    meffect!(state, msr_prob) = measuresitesrandomly!(state, msrop, msr_prob; ITensors_apply_kwargs...)
    return meffect!
end