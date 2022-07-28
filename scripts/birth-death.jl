using DrWatson, Revise
@quickactivate

includet(srcdir("non-parametric.jl"))

using BridgeEstimators
using LaTeXStrings, Plots; gr(dpi = 500)

isgrowing(x; δ = 10) = x[end] > x[1] + δ

λ, μ = .95e-1, 1.e-1
η = [λ, μ]
X₀ = 10
P = BirthDeathProcess(η, X₀)

paths = sample(100, P, 1., 100)
isdangerous = vec(mapslices(isgrowing, paths.values; dims = 1))
dangerous = subsample(paths, isdangerous)
mild = subsample(paths, .!isdangerous)

Rt(paths)

# Sample paths, with dangerous ones highlighted
plot(
    mild.values,
    legend = false,
    color = :gray,
    alpha = .5,
    xlabel = "Observation time",
    ylabel = "Size"
)
plot!(
    dangerous.values,
    legend = false,
    color = :orange,
    linewidth = 2
)
savefig(plotsdir("BD-paths"))

plot(Rt(mild; window_length = 10); legend = false, color = :gray, alpha = .5, label = false, ylims = (.6, 1.3))
plot!(Rt(dangerous; window_length = 10); legend = false, color = :orange, label = false, linewidth = 2)
hline!([λ/μ], label = "True " * L"R_0", color = :blue, linewidth = 2, legend = :bottomright)
savefig(plotsdir("BD-Rt"))


# Compared likelihood functions for process and bridge
steps = 0.1e-1:4e-3:2e-1
etas = [[x, y] for x in steps, y in steps]

for i in 1:sum(isdangerous)
    X = subsample(dangerous, i)
    p1 = heatmap(
        steps, 
        steps,
        [likelihood(BirthDeathProcess(η, 0), X; bridge = false, log = false) for η in etas]',
        legend = false,
        xlabel = "Growth rate λ",
        ylabel = "Death rate μ",
        title = "Likelihood under original process"
    )
    scatter!([η[1]], [η[2]], label = false)
    
    p2 = heatmap(
        steps, 
        steps,
        [likelihood(BirthDeathProcess(η, 0), X; bridge = true, log = false) for η in etas]',
        legend = false,
        xlabel = "Growth rate λ",
        ylabel = "Death rate μ",
        title = "Likelihood under bridge process"
    )
    scatter!([η[1]], [η[2]], label = false)
    
    plot(p1, p2, size = (1000, 600))
    savefig(plotsdir("BD-likelihood" * string(i)))  
end
