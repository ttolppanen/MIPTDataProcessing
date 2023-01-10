# using QuantumOperators

# s : state; a state from the calculated time-evolution
function get_observables()
    observables = Dict(
        :entanglement_half => (; exact = "entanglement(d, L, s, Int(floor(L/2)))", mps = "entanglement(s, Int(floor(L/2)))")
    )
    return observables
end