# using QuantumOperators

# sp : simulation parameters; these include d, L, dt etc.
# s : state; a state from the calculated time-evolution
function get_observables()
    observables = Dict(
        :entanglement_half => (; exact = (s; sp...) ->  entanglement(sp[:d], sp[:L], s, Int(floor(sp[:L]/2))), mps = (s; sp...) ->  entanglement(s, Int(floor(sp[:L]/2))))
    )
    return observables
end