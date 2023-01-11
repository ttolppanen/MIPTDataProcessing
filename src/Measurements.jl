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
    meffect!(state) = measuresitesrandomly!(state, msrop, sp[:p])
    return meffect!
end
function std_n_mps(ITensors_apply_kwargs; sp...)
    msrop = measurementoperators(nop(sp[:d]), siteinds(sp[:state0]))
    meffect!(state) = measuresitesrandomly!(state, msrop, sp[:p]; ITensors_apply_kwargs...)
    return meffect!
end