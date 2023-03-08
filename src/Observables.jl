# using QuantumOperators
# using ITensors
# using LinearAlgebra

ITensors.op(::OpName"n_squared_MIPT_DP", ::SiteType"Boson", d::Int) = Matrix(nop(d)^2)

# sp : simulation parameters; these include d, L, dt etc.
# s : state; a state from the calculated time-evolution
function get_observables()
    observables = Dict(
        :entanglement_half => (; exact = entanglement_half, mps = entanglement_half_mps),
        :mutual_information_first_to_half => (; exact = mutual_information_first_to_half, NOT_IMPLEMENTED_YET = (s; sp...) ->  "NOT_IMPLEMENTED_YET"),
        :mutual_information_half_centered => (; exact = mutual_information_half_centered, NOT_IMPLEMENTED_YET = (s; sp...) ->  "NOT_IMPLEMENTED_YET"),
        :mutual_information_half_centered_region => (; exact = mutual_information_half_centered_region, NOT_IMPLEMENTED_YET = (s; sp...) ->  "NOT_IMPLEMENTED_YET"),
        :half_n_fluctuations => (; exact = half_n_fluctuations, mps = half_n_fluctuations_mps),
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

function mutual_information_half_centered(; sp...)
    a = ceil(Int, sp[:L] / 2) - floor(Int, sp[:L] / 4) # ceil is for the odd-case
    b = sp[:L] - (a - 1)
    return s -> entanglement(sp[:d], sp[:L], s, [a]) + entanglement(sp[:d], sp[:L], s, [b]) - entanglement(sp[:d], sp[:L], s, [a, b])
end
function mutual_information_half_centered_mps(; sp...)
    half = Int(floor(sp[:L]/2))
    return s -> entanglement(s, [1]) + entanglement(sp[:d], sp[:L], s, [half + 1]) - entanglement(sp[:d], sp[:L], s, [1, half + 1])
end

function mutual_information_half_centered_region(; sp...)
    a = ceil(Int, sp[:L] / 2) - floor(Int, sp[:L] / 4) # ceil is for the odd-case
    b = sp[:L] - (a - 1)
    region_length = floor(sp[:L] / 8)
    a = a:(a + region_length)
    b = (b - region_length):b
    return s -> entanglement(sp[:d], sp[:L], s, a) + entanglement(sp[:d], sp[:L], s, b) - entanglement(sp[:d], sp[:L], s, union(a, b))
end
function mutual_information_half_centered_region_mps(; sp...)
    half = Int(floor(sp[:L]/2))
    return s -> entanglement(s, [1]) + entanglement(sp[:d], sp[:L], s, [half + 1]) - entanglement(sp[:d], sp[:L], s, [1, half + 1])
end

function boson_half(; sp...)
    return s -> bosonmean(sp[:d], sp[:L], s, 1:Int(floor(sp[:L]/2)))
end

function half_n_fluctuations(; sp...)
    #NOT_IMPLEMENTED_YET
end
function half_n_fluctuations_mps(; sp...)
    half = Int(floor(sp[:L]/2))
    return s -> sum(expect(s, "n_squared_MIPT_DP"; sites = 1:half)) - sum(expect(s, "N"; sites = 1:half))^2
end