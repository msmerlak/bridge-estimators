using SimpleBirthDeathProcess
using Optim

function loglikelihood(η, X; bridge = false)

    l = loglik(η, X)

    if bridge
        l -= trans_prob(X.state[1], X.state[end], X.k * X.u, η)
    end

    return l
end

function MLE(X; bridge=false, verbose=false, tol=1e-2)

    F(η) = -loglikelihood(η, X; bridge=bridge)

    lower = [1e-5, 1e-5]
    upper = [Inf, Inf]
    initial = ones(2)

    try
        sol = optimize(F, lower, upper, initial, Fminbox(Optim.LBFGS()))
        return sol.minimizer

    catch e
        return missing
    end
end
