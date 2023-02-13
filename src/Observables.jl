# using QuantumOperators

# sp : simulation parameters; these include d, L, dt etc.
# s : state; a state from the calculated time-evolution
function get_observables()
    observables = Dict(
        :entanglement_half => (; exact = entanglement_half, mps = entanglement_half_mps),
        :mutual_information_first_to_half => (; exact = mutual_information_first_to_half, NOT_IMPLEMENTED_YET = (s; sp...) ->  "NOT_IMPLEMENTED_YET"),
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

function mutual_information_first_to_half(; sp...)
    half = Int(floor(sp[:L]/2))
    return s -> entanglement(sp[:d], sp[:L], s, [1]) + entanglement(sp[:d], sp[:L], s, [half + 1]) - entanglement(sp[:d], sp[:L], s, [1, half + 1])
end
function mutual_information_first_to_half_mps(; sp...)
    half = Int(floor(sp[:L]/2))
    return s -> entanglement(s, [1]) + entanglement(sp[:d], sp[:L], s, [half + 1]) - entanglement(sp[:d], sp[:L], s, [1, half + 1])
end

function boson_half(; sp...)
    return s -> bosonmean(sp[:d], sp[:L], s, 1:Int(floor(sp[:L]/2)))
end