# Computing-with-Confidence

Fault Tree Analysis with Imprecise Boolean Logic and Probability Bounds Analysis

`ProbabilityBoundsAnalysis.jl` is a registered Julia package designed by [Dr Ander Gray](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl), and `UncLogic.jl` is a registered Julia package co-designed by [Enrique Miralles-Dolz and Dr. Ander Gray](https://github.com/Institute-for-Risk-and-Uncertainty/UncLogic.jl)
The latest releases can be installed via the Julia package manager:

For the `ProbabilityBoundsAnalysis.jl` package:
```julia
julia> ]
(v1.0) pkg> add ProbabilityBoundsAnalysis
julia> using ProbabilityBoundsAnalysis
```

For the `UncLogic.jl` package:
```
> git clone https://github.com/AnderGray/UncLogic.jl.git
> cd UncLogic
> julia --project
```
```julia
julia> using UncLogic
```

### Related packages:
* [pba.r](https://github.com/ScottFerson/pba.r): R version of this software.
* [RAMAS® RiskCalc](https://www.ramas.com/riskcalc): a commerical software for distribution-free risk analysis.
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): the interval arithmetic package used in this software.
