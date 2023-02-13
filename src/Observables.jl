# using QuantumOperators

# sp : simulation parameters; these include d, L, dt etc.
# s : state; a state from the calculated time-evolution
function get_observables()
    observables = Dict(
        :entanglement_half => (; exact = entanglement_half, mps = entanglement_half_mps),
        :boson_half => (; exact = boson_half, NOT_IMPLEMENTED_YET = (s; sp...) ->  "NOT_IMPLEMENTED_YET")
    )
    return observables
end

function entanglement_half(; sp...)
    return s -> entanglement(sp[:d], sp[:L], s, Int(floor(sp[:L]/2)))
end
function entanglement_half_mps(; sp...)
    return s -> entanglement(s, Int(floor(sp[:L]/2)))
end

function boson_half(; sp...)
    return s -> bosonmean(sp[:d], sp[:L], s, 1:Int(floor(sp[:L]/2)))
end