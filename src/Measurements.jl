# using QuantumOperators
# using QuantumTimeEvolution

function get_measurements()
    measurements = Dict(
        :std_n => (; exact = "std_n_exact(d, L, msr_prob)", mps = "std_n_mps(state0, d, msr_prob)")
    )
    return measurements
end

function std_n_exact(; sp...)
    msrop = measurementoperators(nop(sp[:d]), sp[:L])
    meffect!(state) = measuresitesrandomly!(state, msrop, sp[:msr_prob])
end
function std_n_mps(; sp..., state0, d, msr_prob; ITensors_apply_kwargs...)
    ITensors_apply_kwargs = filter(kv -> kv.first ∉ [:d, :state0, :msr_prob] , sp)
    msrop = measurementoperators(nop(sp[:d]), siteinds(sp[:state0]))
    meffect!(state) = measuresitesrandomly!(state, msrop, sp[:msr_prob]; ITensors_apply_kwargs...)
end