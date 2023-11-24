# Computing with Confidence

The objective of the research work is to demonstrate the ability to propagate confidence structures, in the form of confidence-boxes (c-boxes), with the aim of geenrating generalised confidence bounds over the reliability estimates for engineering systems.
In the work, the relationship between Boolean logic oeprations with system configurations (i.e., series or parallel) is established and the computations no longer just assumes independence, but also the uncertainties in the dependencies between components and path sets (i.e., dual-level dependencies).
To perform the necessary computations, the Julia progamming is used with the following packages: `ProbabilityBoundsAnalysis.jl` and `UncLogic.jl`. The role of the repository is to provide a platform to present the codes which serve as a tutorial for those who are interested in learning and implementing the proposed methods presented in the literature.

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
