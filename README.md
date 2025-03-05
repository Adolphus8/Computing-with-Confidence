# Computing with Confidence

The objective of the research work is to demonstrate and evaluate the feasibility to propagate confidence structures, in the form of non-consonant confidence-boxes (c-boxes), with the aim of generating generalised confidence bounds over the reliability estimates for engineering systems.
In the work, the relationship between Boolean logic oeprations with system configurations (i.e., series or parallel) is established and the computations no longer just assumes independence, but also the uncertainties in the dependencies between components.
To perform the necessary computations, the R progamming code is used with the following packages: `pba BETTER.r`. The role of the repository is to provide a platform to present the codes which serve as a tutorial for those who are interested in learning and implementing the proposed methods presented in the literature.

`pba BETTER.r` is a R programming package to for probability bounds analysis while the `sra.r` is an educational R programming package for risk assessment, both designed by [Prof. Scott Ferson](https://github.com/ScottFerson/pba.r.git).

## 1) Coherent systems:

The repository consists of the codes to the Illustrative example as well as the three case studies presented in the paper. The following codes are:
1) Illustrative Examples.r
2) Case study 1 - Pressure tank system.r
3) Case study 2 - TRIGA Research Reactor.r
4) Case study 3 - Bridge structure.r
5) Images.m (Run this MATLAB code last as it is a file to generate the images presented in the manuscript)
6) Supplementary folder - Containing MATLAB codes to the three case studies which seeks to compare the Structure function solution to that of the Boolean expression using Monte Carlo (under independence assumption only)

### Related packages:
* [ProbabilityBoundsAnalysis.jl](https://github.com/AnderGray/ProbabilityBoundsAnalysis.jl): Julia version of this software.
* [RAMAS® RiskCalc](https://www.ramas.com/riskcalc): a commerical software for distribution-free risk analysis.
* [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl): the interval arithmetic package used in this software.

## Reference:
* A. Lye, W. Vechgama, M. Sallak, S. Destercke, S. Ferson, and S. Xiao (2024). Advances in the reliability analysis of coherent systems under limited data with Confidence boxes. *ASCE-ASME Journal of Risk and Uncertainty in Engineering Systems Part A: Civil Engineering* (Accepted). doi: TBC
