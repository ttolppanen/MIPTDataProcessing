# using QuantumStates

export get_initial_states

function get_initial_states()
    states = Dict(
        :zeroone => (; exact = (; sp...) -> zeroone(sp[:d], sp[:L]), mps = (; sp...) -> zeroonemps(sp[:d], sp[:L])),
        :onezero => (; exact = (; sp...) -> onezero(sp[:d], sp[:L]), mps = (; sp...) -> onezeromps(sp[:d], sp[:L]))
    )
    return states
end
