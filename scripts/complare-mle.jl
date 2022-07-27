using DrWatson, Revise
@quickactivate
includet(srcdir("mle-birth-death.jl"))

using Plots; gr(dpi = 500)

isgrowing(X; δ = 10) = X.state[end] > X.state[1] + δ


η = [.95e-1, 1.e-1]
paths = [rand_discrete(5, 50, 1, η) for _ in 1:100]
dangerous = paths[isgrowing.(paths)]
mild = paths[.!isgrowing.(paths)]


# Sample paths, with dangerous ones highlighted
plot(
    [g.state for g in mild],
    legend = false,
    color = :gray,
    alpha = .5,
    xlabel = "Observation time",
    ylabel = "Size"
)
plot!(
    [g.state for g in dangerous],
    legend = false,
    color = :orange,
    linewidth = 2
)
savefig(plotsdir("BD-paths"))


# Compared likelihood functions for process and bridge
steps = 0.1e-1:4e-3:2e-1
etas = [[x, y] for x in steps, y in steps]

for i in 1:length(dangerous)
    X = dangerous[i]
    p1 = heatmap(
        steps, 
        steps,
        [exp(loglikelihood(η, X; bridge = false)) for η in etas]',
        legend = false,
        xlabel = "Growth rate λ",
        ylabel = "Death rate μ",
        title = "Likelihood under original process"
    )
    scatter!([η[1]], [η[2]], label = false)
    
    p2 = heatmap(
        steps, 
        steps,
        [exp(loglikelihood(η, X; bridge = true)) for η in etas]',
        legend = false,
        xlabel = "Growth rate λ",
        ylabel = "Death rate μ",
        title = "Likelihood under bridge process"
    )
    scatter!([η[1]], [η[2]], label = false)
    
    plot(p1, p2, size = (1000, 600))
    savefig(plotsdir("BD-likelihood" * string(i)))  
end
    

