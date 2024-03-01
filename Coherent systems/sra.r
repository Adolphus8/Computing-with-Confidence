##########################################################################
# Monte Carlo Simulation S4 Library for the R Language                
#                                                                          
# Institute for Risk and Uncertainty, University of Liverpool, 
# Copyright 2023 Scott Ferson; All rights reserved worldwide     
# ferson(at)liverpool(dot)ac(dot)uk, or sandp8(at)gmail(dot)com
##########################################################################
#
# Place this file on your computer and, from within R, select 
# File/Source R code... from the main menu.  Select this file and
# click the Open button to read it into R.  You should see the 
# completion message ":sra> library loaded".  Once the library 
# has been loaded, you can define probability distributions 
#
#       a = normal(5,1)
#       b = uniform(2,3)
#
# and perform mathematical operations on them, including
# convolutions such as 
#
#       a  +  b
#
# or generalized convolutions such as 
#
#       a + b * beta(2,3) / a 
#
# Many different distribution shapes are supported via the following
# functions
#
#       bernoulli, beta (B), binomial (Bin), cauchy, chi, chisquared, 
#       delta, dirac, discreteuniform, exponential, exponentialpower, 
#       extremevalue, F, f, fishersnedecor, fishertippett, fisk, frechet, 
#       gamma, gaussian, geometric, generalizedextremevalue (GEV), 
#       generalizedpareto, gumbel, histogram, inversegamma, laplace, 
#       logistic, loglogistic, lognormal (L), logtriangular, loguniform, 
#       negativebinomial, normal (N), pareto, pascal, powerfunction, 
#       poisson, quantiles, rayleigh, reciprocal, shiftedloglogistic, 
#       skewnormal (SN), student, trapezoidal, triangular (T), 
#       uniform (U), weibull,
#
# as well as several functions for fitting distributions to data using 
# the method of matching moments such as
#
#       MMbernoulli, MMbeta, MMbetabinomial, MMbinomial          
#       MMchisquared, MMdoubleexponential, MMexponential, 
#       MMextremevalue, MMF, MMgamma, MMgaussian, 
#       MMgeometric, MMgumbel, MMlaplace, MMlogistic, 
#       MMlognormal, MMloguniform, MMnormal, MMpareto, 
#       MMpascal, MMpoisson, MMpowerfunction, MMrectangular, 
#       MMstudent, MMt, MMtriangular, MMuniform, MMuniform1, 
#
# or maximum likelihood such as 
#
#       MLbernoulli, MLdoubleexponential, MLexponential, 
#       MLgamma, MLgaussian, MLgeometric, MLlaplace, 
#       MLlognormal, MLnormal, MLpareto, MLpascal, 
#       MLpoisson, MLrectangular, MLuniform, MLweibull, 
#
# or maximum entropy such as 
#
#       MEdiscretemean, MEdiscreteminmax, MEmeansd, 
#       MEmeanvar, MEminmax, MEminmaxmean, 
#       MEminmaxmeansd, MEminmaxmeanvar, MEminmean, 
#       MEmmms, MEquantiles, 
#
# as well as functions useful in Bayesian inference such as
#
#       km, KN, betabinomial, gammaexponential, BCbernoulli, 
#       BCbinomial, BCexponential, BCgeometric, BCuniform, 
#       BCnegativebinomial, BCnormal, BCnormal.knownmu, 
#       BCnormal.knownsigma, BCpareto, BCpoisson.
#
# Note that the quantiles and histogram functions allow you to 
# create distributions of arbitrary shape.
#
# You can use regular mathematical operations such as +, -, *, /, 
# and ^ to compute with these distributions. Comparison operators
# such as <, >, <=, >=, and even == are defined for distributions.
# Comparison operations typically yield probabilities (that is, 
# numerical values between zero and one), rather than the Boolean
# constants FALSE and TRUE.  These probabilities represent the 
# proportion of times the comparions would be true for respective
# values from the distributions.
#
# Using such probabilities, you can do probabilistic logic with 
# the functions
#
#	not, and, or, xor, cond, nand, nor, imply
#
# For instance, the and() function computes the probability of the 
# conjunction or intersection of two events of known probabilities. 
# Likewise, the or() function gives the analogous probability of 
# disjunction or union of two events given their respective 
# probabilities.  Other than the not() function which gives the 
# negation of a single probability, each of these functions expects
# two probabilites and supports an optional third parameter r that 
# specifies the stochastic dependence to be assumed about the input 
# probabilities.  If this parameter is given as zero or is absent, 
# the functions assume stochastic independence.  For instance, 
# and(0.2, 0.3) yields the probability 0.06, while or(0.2, 0.3) 
# gives the so-called probabilistic sum 0.44.  Specify r as +1 to 
# indicate that the dependence should be assumed to be perfect, or 
# -1 for opposite dependence.  The probabilities may themselves be 
# distributions, as in, for instance, the calculation the expression
# and(U(.2,.3), N(.5,.1)).  Note that, in this case, the r parameter 
# for these logical functions will interact with whatever functional 
# dependence there may be between distributional inputs.  Impossible
# values for r are automatically reduced to perfect dependence, or
# increased to opposite dependence, as appropriate.
#
# Several standard and new mathematical transformations and 
# binary functions have also been introduced or extended to handle 
# probability distributions, including
#
#       exp, log, sqrt, abs, round, trunc, ceiling, floor, sign,  sin, cos, tan,
#       asin, acos, atan, atan2, reciprocate, negate, not, complement, 
#       rescale, truncate, pmin, pmax, and, or, xor, nand, nor, cond, imply.
#
# The stochastic dependence between the results of these functions 
# and their input is automatically preserved.  For instance, entering 
# the expressions
#
#	A = N(5,1)
#	cor(A, log(A))
#
# will yield a value around 0.99 as the correlation between the random 
# variable A and its logarithms, even though entering cor(A,log(N(5,1))) 
# will give a value very close to zero.
#
# Likewise, the convolution distributions obtained by combining input 
# distributions with the arithmetic operators +, -, *, /, ^, pmin, pmax 
# will automatically encode the dependencies with the respective inputs.
# For instance, if A = N(5,1), the value returned by cor(A, A+U(10,15)) 
# will be around 0.57.  The correlation increases with the width of the 
# uniform addend decreases.
#
# As explained above, the outputs of the comparison operators <, etc., 
# are reduced to probabilities.  If you wish to capture the dependencies 
# implicit in the underlying distribution of truth values, you can use a
# construction like mc((A@x<B@x)+0), where A and B are the inputs 
# in A<B, which yields a Bernoulli distribution of truth values (with FALSE 
# represented by zero, and TRUE represented by one) having a mean equal 
# to the probability A<B but with the proper dependencies with the random 
# variables A and B.
#
# Separately constructed distributions are assumed to be independent by 
# default.  Distributions from separate calculations are typically also 
# independent of one another unless the calculations share inputs or 
# correlations with other variables.  Function results will be dependent 
# on the variables used to create them.  In particular, assigning a 
# variable containing a probability distribution to another variable makes 
# the two variables perfectly dependent.  To create a copy of the values 
# distribution that is independent of it, you can use the samedistribution 
# function, e.g., C = samedistribution(A). 
#
# You can account for perfect or opposite dependence between distributions 
# by giving their dependence when you construct them with expressions like 
#
#       d = beta(2,5, r=a)  # perfect dependence (comonotonicity)
# or
#       d = beta(2,5, r=-a) # opposite dependence (countermonotonicity)
#
# In the first example, the variable d will be perfectly dependent with 
# the variable a.  In the second case, because of the minus sign in front 
# of the a, the variables d and a will be oppositely dependent.
#
# Perfect and opposite dependencies are automatically mutual, so it is not 
# necessary to explicitly make any reciprocal assignments.  Thus
#
#       a = N(5,1)
#       b = U(2,3, r=a)
#       c = N(15,2, r=b)
#
# suffices to link c with a and vice versa.  The assignments automatically 
# make a, b, and c mutually perfectly dependent.  Naturally, it is not 
# possible to be perfectly (or oppositely) dependent on more than one 
# quantity unless they are also mutually dependent in the same way.  
#
# It is also possible to specify intermediate correlations in between 
# opposite and perfect (that are not independent) although the methods for 
# doing this are approximate, so you should always check that your realized 
# correlations are close to the correlations you planned.  To introduce
# correlation to a pair of distributions, first construct one distribution, 
# and then use the r option when constructing the second distribution to 
# specify the planned correlation.  The value for r should be a pair made
# with the c() function consisting of a (Pearson product-moment) correlation 
# coefficient, which is some value between -1 and +1, and the name of the 
# (first) variable with which the new distribution should be correlated. 
# For example,
#
#	a = N(5,1)
#	b = beta(2,18, r=c(0.8,a))
#	cor(a,b)   
#	plot(a,b)
#
# creates correlated random variables, revealed by the plot() function to 
# be far from independent.  In this case, the observed correlation will be 
# around 0.77.  (As already mentioned, the software uses approximate methods 
# to model correlations.  Increasing the correlation coefficient to 0.83 in 
# this case produces a correlation closer to 0.8.  Values outside the range
# [-1,+1] generate an error.)
#
# When specifying correlations pairwise is cumbersome or insufficient, you
# call also construct multiple random variables with a given matrix of cross 
# correlations.  You will need to installing the MASS package for R.  Simply
# specify intermediate dependencies among several variables in a matrix 
# of Pearson correlation coefficients, and call the correlations() function 
# to produce "r" vectors for use in constructing each variable, deploying 
# them as in the following example:
#
#	# planned correlation matrix
#	s = c(  1.00,  0.60, -0.30,
#	           0.60,  1.00, -0.40,
#	          -0.30, -0.40,  1.00)
# 
#	correlations(s)   # read corr matrix s to generate "r" vectors
# 
#	w = uniform(2,5, r=MC$r[,1])    
#	x = poisson(5,   r=MC$r[,2])
#	y = gumbel(2,3,  r=MC$r[,3])
# 
#	# pairwise bivariate plots
#	plotcorrs(c(w,x,y),names=c('w','x','y'))  
#	corrs(c(w,x,y))
# 
# Dependence between variables can have dramatic effects on any calculations 
# that involve these variables.  Typically the Monte Carlo simulations take 
# account of the dependence automatically.  You can alternatively force the 
# software to assume two variables to have perfect or opposite dependence when 
# you do operations on them.  For instance, you may assume perfect dependence 
# with an expression like
#
#       perfectconv.mc(a, b, '*')
#
# or assume oposite dependence with an expression like
#
#       oppositeconv.mc(a, b, '+')
#
# After making the calculations, however, a and b are not correlated with each 
# other, nor will they natural correlation with the output distribution.
#
# For convenience, and perhaps because life is too drab, many new commands for 
# plotting distributions are defined, such as
# 
#       edf, black, blue, brick, brown, chartreuse, cyan, gray, green,     
#       khaki, magenta, navy, olivedrab, orange, pink, purple, red, 
#       rust, salmon, sienna, tan, teal, white, yellow,
#
# as well as their three-letter abbreviations. Also, several standard commands 
# allow you to see the resulting distributions, such as
#
#       show(c)
#       summary(d)
#       plot(a)
#       lines(b, col='blue')
#
# By default, simply entering a distribution's name on the R Console will 
# display it graphically. Set MC$plotting to FALSE to disable this.
#
# There are a variety of standard and new functions you can use with
# distributions to characterize them, such as 
#
#       mean, sd, var, median, quantile, fivenum, left, right, prob, cut,
#       percentile, iqr, IQR, random, randomrange, specific, psd, pvar,
#       cv, skewness, kurtosis, ci, qci, KSbands,
#
# and several functions to compare distributions, including
#
#	areametric, smirnov, cor, dom.
#
# Several mass reassignment functions are supported such as 
#
# 	below, above, lowest, highest, smin, smax, between,
#
# and miscellaneous constructors such as 
#
# 	mixture, selfsum, extremise.
#
# Note that, unlike transformations like exp and log, these functions 
# refer to distributions and not the variates themselves. Consequently, 
# they do not preserve any implicit dependencies, so selfsum(a,2) is 
# the same as a+shuffle(a), which leads to a distribution with a smaller 
# variance than that of a+a = a*2.
#
# Programming in R is imperfect.  Although distributions can be mixed 
# together with real values in many operations, this is not possible with 
# some functions including the two-argument log.  If you need to use 
# real values among the arguments for these functions, it may help to
# wrap the real values in the mc function. Although log(5, normal(5,1)) 
# precipitates an error, the expression log(mc(5), N(5,1)) will work as 
# expected.
# 
##########################################################################
# 
# END OF NOTES FOR USER 
#
# More notes for programmers are in comments throughout the code below
#
##########################################################################


##########################################################################
# Global constants #
##########################################################################
MC <- new.env()
MC$many <- 20000                # Monte Carlo replications (increase for better precision)
MC$lwd <- 3                           # line width for plotting
MC$plotting <- TRUE              # if TRUE, distributions are plotted whenever they are "show"n
MC$plotting.every <- FALSE    # if TRUE, every distribution that's created is automatically plotted
MC$cumulative <- TRUE          # if TRUE, plot CDFs, if FALSE, plot CCDFs (exceedance risks)
MC$distribs <- 0                      # number of Monte Carlo distributions that have been created
MC$verbose <- 2                     # how much warning messaging is wanted
MC$r <- 0                               # place for the correlations function to put its correlated sequences
RStudio <- Sys.getenv("RSTUDIO") == '1'

# minimize footprint
if ("quiet" %in% ls()) MC$scott.used.to.be.quiet <- quiet

####################################
# Some ancillary utility functions 
####################################

Gamma <- gammafunc <- function(x) .Primitive("gamma")(x)  # 'gamma' is redefined to be a distribution

ami <- function(i,j) isTRUE(all.equal(i,j)) 

##########################################################################
# Class specification and general MC distribution constructor 
##########################################################################

quiet <- setClass('mc',representation(x='numeric', n='integer', bob='integer', id='character')) 

uniquemc <- function() { MC$distribs <- MC$distribs + 1; return(paste(MC$distribs)) }

mc <- function(x, perfect=NULL, opposite=NULL, bob=NULL, interpolation='linear') {
  if (missing(x)) x <- rep(0,MC$many) 
  if (inherits(x,'mc')) p <- x  else
    {
      if (!is.number(x)) stop('Monte Carlo distribution must be numerical') 
      if (length(x) != MC$many) x <- interpolate(x,interpolation); 
      unique <- uniquemc()
      id <- paste('MC',unique,sep='')
      if (!missing(bob)) unique <- bob
      p <- new("mc", x=x, n=as.integer(MC$many), bob=as.integer(unique), id=id)
    }
  class(p) <- c('mc')
  if (MC$plotting.every) try(plot.mc(p))
  p
  }

as.mc <- as.mc.numeric <- as.mc.interval <- mc

makemc <- function(...) {  #  works for a bunch of args separated by commas, or for a single list of args
  elts <- list(...)
  if (mode(elts[[1]])=='list') elts <- elts[[1]]
  if (ami(1,length(elts))) return(mc(...))
  for (i in seq(along=elts)) if (!is.mc(elts[[i]])) elts[[i]] <- mc(elts[[i]])
  elts
  }


##########################################################################
# Interpolation schemes            
##########################################################################

# The mc() constructor accepts lists of x-values to define
# the distribution.  Both linear and spline interpolation 
# schemes can be used with these lists. The spline 
# interpolation scheme requires the splines library to be 
# loaded and needs at least 4 values.  
#    See also the quantiles() constructor when you have
# quantiles (or percentiles or fractiles).

interpolate <- function(u, interpolation='linear') switch(interpolation, "spline"={spline.interpolate(u)}, {linear.interpolate(u)})

# needs at least 4 points and the interpSpline function from the library(splines)
spline.interpolate <- function(u) {library(splines); return(predict(interpSpline(1:length(u),u),nseg=MC$many-1)$y)}  # needs the splines library

linear.interpolate <- function(V) {
  m <- length(V) - 1
  if (ami(m,0)) return(rep(V,MC$many))
  if (ami(MC$many,1)) return(c(min(V),max(V)))
  d <- 1/m
  n <- round(d * MC$many * 20)
  if (n==0) c = V else {
  c <- NULL
  for (i in 1:m) {
    v <- V[[i]]     
    w <- V[[i+1]]   
    c <- c(c, seq(from=v, to=w, length.out=n))
    } }
  u <- rep(0,MC$many)
  for (k in 1:MC$many) u[[k]] <-c[[1+round((length(c)-1)*(k-1)/(MC$many-1))]]
  u
  }

##########################################################################
# Access methods #
##########################################################################

left       <- function(x,...) {UseMethod("left")};         left.default    <- function(x) return(min(x))
right      <- function(x,...) {UseMethod("right")};        right.default   <- function(x) return(max(x))
min        <- function(...,na.rm=FALSE){UseMethod("min")}; min.default     <- function(..., na.rm = FALSE) base::min(..., na.rm = na.rm)
max        <- function(...,na.rm=FALSE){UseMethod("max")}; max.default     <- function(..., na.rm = FALSE) base::max(..., na.rm = na.rm)
#pmin       <- function(...,na.rm=FALSE){UseMethod("pmin")};pmin.default    <- function(..., na.rm = FALSE) base::pmin(..., na.rm = na.rm)
#pmax       <- function(...,na.rm=FALSE){UseMethod("pmax")};pmax.default    <- function(..., na.rm = FALSE) base::pmax(..., na.rm = na.rm)
pmin <- function(...,na.rm=FALSE) {UseMethod("pmin")};     pmin.default <- function(..., na.rm = FALSE) base::pmin(...,na.rm=na.rm)
pmax <- function(...,na.rm=FALSE) {UseMethod("pmax")}; pmax.default <- function(..., na.rm = FALSE) base::pmax(...,na.rm=na.rm)
ci       <- function(x...,c=0.95,alpha=(1-c)/2,beta=1-(1-c)/2,s = sort(x),na.rm = FALSE){UseMethod("ci")};ci.default    <- function(x,..., c=0.95, alpha=(1-c)/2, beta=1-(1-c)/2, s = sort(x), na.rm = FALSE)  c(s[round(alpha*length(x))],s[round(beta*length(x))])
prob       <- function(x,s,...) {UseMethod("prob")}
alpha      <- function(x,s,...) {UseMethod("prob")}
shuffle    <- function(x) {UseMethod("shuffle")}; shuffle.default <- function(x) x[order(runif(length(x)))]
iqr          <- function(x) {UseMethod("iqr")}
IQR        <- function(x, na.rm = FALSE, type = 7) {UseMethod("IQR")}; IQR.default <- function(x, na.rm = FALSE, type = 7) stats::IQR(x,na.rm,type)
fivenum   <- function(x, na.rm = TRUE) {UseMethod("fivenum")}; fivenum.default <- function(x, na.rm = TRUE) stats::fivenum(x,na.rm)
shift        <- function(ss,x) {UseMethod("shift")};shift.default <- function(ss,x) ss*x
mult       <- function(m,x) {UseMethod("mult")};mult.default <- function(m,x) m*x
atan2      <- function(y,x) {UseMethod("atan2")};atan2.default    <- function(y,x) base::atan2(y,x)
fractile     <- function(x,p,...) {UseMethod("fractile")}
percentile <- function(x,p,...) {UseMethod("percentile")}
negate     <- function(x,...) {UseMethod("negate")};       negate.default  <- function(x,...) return(-x)
reciprocate<- function(x,...) {UseMethod("reciprocate")};  reciprocate.default  <- function(x,...) return(1/x)
complement <- function(x,...) {UseMethod("complement")};   complement.default  <- function(x,...) return(1-x)
#and        <- function(x,y,...) {UseMethod("and")};        and.default  <- function(x,y,...) and.mc(x,y)
#or         <- function(x,y,...) {UseMethod("or")};         or.default   <- function(x,y,...) or.mc(x,y)
and        <- function(x,y,r=0,...) {UseMethod("and")};        and.default  <- function(x,y,r=0,...) and.numeric(x,y,r)
or         <- function(x,y,r=0,...) {UseMethod("or")};         or.default   <- function(x,y,r=0,...) or.numeric(x,y,r)
not        <- function(x,...)   {UseMethod("not")};        not.default  <- function(x,...)   complement.mc(x)
cor          <- function(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman")) {UseMethod("cor")}; cor.default <- function(x, y = NULL, use = "everything", method = c("pearson", "kendall", "spearman")) stats::cor(x,y,use,method)
#cut         <- function(x,s,...) {UseMethod("fractile")}
#trunc      <- function(x,m,M,...) {UseMethod("truncate")}; truncate.default<- function(x,...) return(base::truncate(x,...))
#round      <- function(x,...) {UseMethod("round")};        round.default   <- function(x,...) return(base::round(x,...))
#ceiling    <- function(x,...) {UseMethod("ceiling")};      ceiling.default <- function(x,...) return(.Primitive("ceiling")(x,...))
#sign       <- function(x,...) {UseMethod("sign")};         sign.default    <- function(x,...) return(.Primitive("sign")(x,...))

reps <- function(x) if (class(x)=='mc') if (x@n==length(x@x)) x@n else "Internal step count error" else length(x)

min.mc <- left.mc <- function(x) base::min(x@x)

max.mc <- right.mc <- function(x) base::max(x@x)

# function revised to be unary transformation
#sign.mc <- function(x) if (right(x)<0) -1 else if (0<left(x)) +1 else if (is.scalar(x) && identical(x@x[[1]],0)) 0 else if (right(x)<=0) c(-1,0) else if (0<=left(x)) c(0,1) else c(-1,1)

# range.mc <- function(x, na.rm = FALSE) base::range(x@x)  # hard to believe I need to have the na.rm argument!
range.mc = function(..., na.rm = FALSE) {  # works for a bunch of args separated by commas, or for a single list of args
  elts <- list(...)
  if (mode(elts[[1]])=='list') elts <- elts[[1]]
  X = NULL
  for (e in elts) if (is.mc(e)) X = c(X, base::range(e@x)) else X = c(X, base::range(e))
  return(range(X))
  }
#a = N(5,7)
#b = a^2
#range(1,3)  # 1 3
#range(c(4,7))  # 4 7
#range(a)  # -26.2736  30.2651
#range(a,b)  #  -26.2736 915.9760
#range(c(a,b))  # Error in min(x, na.rm = na.rm) : invalid 'type' (list) of argument
#range.mc(c(a,b))  # -26.2736 915.9760
#rbyc(5,1)
#pl(1,3)
#pl(c(4,7))
#pl(a)
#pl(a,b)
#pl(c(a,b))

ci.mc <- function(x, c=0.95, alpha=(1-c)/2, beta=1-(1-c)/2, s = sort(x@x)) c(s[round(alpha*x@n)],s[round(beta*x@n)])

iqr.mc <- function(x) c(percentile.mc(x,0.25), percentile.mc(x,0.75))

IQR.mc <- function(x) diff(iqr.mc(x))

cut.mc <- quantile.mc <- fractile.mc <- percentile.mc <- function(x, p, tight=TRUE) {
  if (any((p<0) | (1<p))) warning('Second argument must be a probability between zero and one')
  sort(x@x)[round(x@n*p)]
  }

#random <- function(x,n=1) x@x[shuffle(1:length(x@x))[1:n]]
#random <- function(x,n=1) if (class(x) == 'numeric') return(rep(x,n)) else return(x@x[shuffle(1:length(x@x))[1:n]])
random <- function(x,n=1) if (class(x)=='mc') return(x@x[sample.int(length(x@x),n,replace=TRUE)]) else return(sample(x,n,replace=TRUE)) 

randomrange <- function(x,n=1) {L = left(x); return(runif(n) * (right(x) - L) + L)}

specific <- function(x) random(x,1)

prob.mc <- alpha.mc <- function(x, s) {
#  x <- makemc(x)
  length(x@x[x@x<s])/length(x@x)
  }

xprob.mc <- function(x, s) {    # required by the inequality comparisons
#  x <- makemc(x); 
  length(x@x[x@x<=s])/length(x@x)
  }
  
density.mc <- function(x, ...) density(x@x, ...)  
  
right.mc <- function(x) max(x@x)

left.mc <- function(x) min(x@x)

#median.mc <- function(x, na.rm = FALSE) percentile.mc(x,0.5)
median.mc <- function(x) percentile.mc(x,0.5) # we're neglecting the na.rm parameters (there shouldn't be any NAs)

mean.mc <- function(x) mean(x@x)             
var.mc <- function(x) var(x@x)
sd.mc <- function(x) sd(x@x)   

psd = function (x) { # population standard deviation, i.e., not sample standard deviation
  n = reps(x)
  sqrt(((n-1)*sd(x)^2)/n)
  }

pvar = function(x) { # population variance, i.e., not sample variance
  n = reps(x)
  ((n-1)*var(x))/n
  }

cv = function(x) return(sd(x)/mean(x))

skewness = function (x, type = 3) {  #  g1, G1, ?, Galton, Pearson 2  # modified from the R library e1071 and the NIST website 
  if (class(x)=='mc') x = x@x
  n = length(x)
  cx = x - mean(x)
  y = sqrt(n) * sum(cx^3)/(sum(cx^2)^(3/2))
  f = fivenum(x);                           
  return(switch(type, y,   y * sqrt(n * (n - 1))/(n - 2), y * ((1 - 1/n))^(3/2), (f[2] + f[4] - 2* f[3])/(f[4] - f[2]), 3 * (mean(x) - f[3])/psd(x)))
  }

kurtosis = function (x, type = 3) {  # type 1 is K1 = k1-3,     type=4 is k1 = Kurt (raw kurtosis, rather than excess kurtosis)
  if (class(x)=='mc') x = x@x
  n = length(x)
  x = x - mean(x)
  r = n * sum(x^4)/(sum(x^2)^2)
  return(switch(type, r-3, ((n + 1) * (r - 3) + 6) * (n - 1)/((n - 2) * (n - 3)), r * (1 - 1/n)^2 - 3, r))
  }

#g1 = function(y) {               # skewness(type=1)      # https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
#  n = length(y)
#  m = mean(y)
#  (sum(((y-m)^3)/n))/(psd(y)^3)
#  }
#G1 = function(y) {              # skewness(type=2)      # https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
#  n = length(y)
#  (sqrt(n*(n-1)) / (n-2)) * g1(y)
#  }
#k1 = function(y) {             # kurtosis(type=4)        # https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
#  n = length(y)
#  m = mean(y)
#  (sum(((y-m)^4)/n))/(psd(y)^4)
#  }
#K1 = function(y) k1(y)-3  # kurtosis(type=1)         # https://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
#
#g2 = function(x) { cx = x - mean(x); n = length(x); (((sum(cx^4))/n)/((sum(cx^2)/n)^2)) -3 } # identical to K1(); https://en.wikipedia.org/wiki/Kurtosis#A_natural_but_biased_estimator
#
#G2 = function(x) { n = reps(x);  ((n-1) / ((n-2)*(n-3))) * ((n+1)*g2(x) + 6) } # https://en.wikipedia.org/wiki/Kurtosis#Standard_unbiased_estimator
#
#x = c(0, 3, 4, 1, 2, 3, 0, 2, 1, 3, 2, 0, 2, 2, 3, 2, 5, 2, 3, 999) # https://en.wikipedia.org/wiki/Kurtosis#A_natural_but_biased_estimator
#skewness(x,type=1);  g1(x)            # 4.129251     # this is Excel's skew.p() function
#skewness(x,type=2);  G1(x)            # 4.471885    # this is Excel's skew() function
#skewness(x,type=3)                        # 3.823462
#skewness(x,type=4)                        # 0.3333333
#skewness(x,type=5)                        # 0.689689
#skewness(x); skewness(x,type=3)   # 3.823462
#kurtosis(x,type=1); K1(x)               # 15.05143 
#kurtosis(x,type=2); G2(x)               # 19.99843   # this is Excel's kurt() function
#kurtosis(x,type=3)                          # 13.29141
#kurtosis(x,type=4); k1(x)               # 18.05143
#kurtosis(x); kurtosis(x)                   # 13.29141
#
#x = rnorm(10,5,1)
#x = c(5.27299941025433, 5.12272209128741, 5.8878308803049, 3.53606966447563, 4.54396345387994, 4.68443023349595, 4.09304301730771, 7.19773966674417, 
#4.06836083160741, 5.34456409541087)
#skewness(x,type=1); g1(x)             # 0.7344248   
#skewness(x,type=2); G1(x)             # 0.8709207   
#skewness(x,type=3)                        # 0.6270629
#skewness(x,type=4)                        # -0.2952769
#skewness(x,type=5)                        # 0.2152019
#skewness(x); skewness(x,type=3)   # 0.6270629
#kurtosis(x,type=1); K1(x)               # 0.09722562
#kurtosis(x,type=2); G2(x)               # 1.136167
#kurtosis(x,type=3)                          # -0.4912473
#kurtosis(x,type=4); k1(x)                # 3.097226 
#kurtosis(x); kurtosis(x,3)                 # -0.4912473
#
#x = rnorm(100000,5,1)
#skewness(x,type=1);  g1(x)            # 0.01387489    
#skewness(x,type=2);  G1(x)            #0.0138751    
#skewness(x,type=3)                        # 0.01387469
#skewness(x,type=4)                        # -0.004071293 
#skewness(x,type=5)                        # -0.0005880126
#skewness(x); skewness(x,type=3)   # 0.01387469
#kurtosis(x,type=1); K1(x)               # -0.02111881
#kurtosis(x,type=2); G2(x)               # -0.02105986
#kurtosis(x,type=3)                          # -0.02117839 
#kurtosis(x,type=4); k1(x)               # 2.978881
#kurtosis(x); kurtosis(x,3)                # -0.02117839

#fivenum.mc <- function(x, na.rm = FALSE) {fn=matrix(c(minimum=left(x), firstquartile=percentile.mc(x,0.25), median=median.mc(x), thirdquartile=percentile.mc(x,0.75), maximum=right(x)),nrow=1); colnames(fn) = c('minimum','1stquartile','median','3rdquartile','maximum'); fn}
fivenum.mc <- function(x, na.rm = FALSE) c(minimum=left(x), firstquartile=percentile.mc(x,0.25), median=median.mc(x), thirdquartile=percentile.mc(x,0.75), maximum=right(x))

#stats.mc <- function(x, na.rm = FALSE) {fn=matrix(c(minimum=left(x), firstquartile=percentile.mc(x,0.25), median=median.mc(x), thirdquartile=percentile.mc(x,0.75), maximum=right(x),mean=mean(x),std=sd(x),var=var(x),sk=skewness(x),kr=kurtosis(x)),nrow=1); colnames(fn) = c('min','1stqtile','median','3rdqtile','max','mean','std','var','skewness','kurtosis'); fn}
stats.mc <- function(x, na.rm = FALSE) c(minimum=left(x), firstquartile=percentile.mc(x,0.25), median=median.mc(x), thirdquartile=percentile.mc(x,0.75), maximum=right(x), mean=mean(x), sd=sd(x), var=var(x), skewness=skewness(x), kurtosis=kurtosis(x), n=reps(x))

cor.mc = function(a,b) stats::cor(a@x,b@x)

qci <- function(a,p=0.5,conf=0.95) { # conf x 100% confidence interval for the pth fractile of the distribution a (see Morgan and Henrion)
 # where Yk  is the (n - k + 1)th largest value from the Monte Carlo simulation,     i = floor(np - b), j = ceiling(np + b), b = zinv((1 - a)/2)sqrt(np(1 - p))
  n = a@n
  b = qnorm((1-conf)/2) * sqrt(n*p*(1 - p))
  i = floor(n*(1-p)-b)
  j = ceiling(n*(1-p)+b)
  s = sort(a@x)
  return(c(s[n-i+1], s[n-j+1]))
  }
## is it reasonable to use qci() with MC distributions, which have arbitrarily large replicate sizes (which don't seem to be the same as statistical "sample sizes")?
#MCSmanySAFE = MC$many
#par(mfrow=c(2,3))
#replace_na = function(x,s) ifelse(is.na(x), s, x)
#for (some in c(50,100,200,500,1000,20000)) {
#  MC$many = some
#  a = N(10,2)
#  plot(a)
#  title(paste(some,'replications'))
#  for (p in 1:99/100) lines(replace_na(qci(a,p),20),rep(p,2),col=2,lwd=2) 
#  }  
#MC$many = MCSmanySAFE

qCI = function(a, p=0.5,conf=0.95) { # Bill Huber's alternative near-symmetric distribution-free confidence interval
  quantile.CI <- function(n, p, alpha=0.05) {  # yields indices into the sorted data for a quantile defined by probability p
    # Notice that the subfunction quantile.CI doesn't depend on the data in the input a, so its outputs could easily be memoised 
    # Searches over a small range of upper and lower order statistics for the closest coverage to 1-alpha (but not less than it, if possible)
    # See https://stats.stackexchange.com/questions/99829/how-to-obtain-a-confidence-interval-for-a-percentile
    u <- qbinom(1-alpha/2, n, p) + (-2:2) + 1
    l <- qbinom(alpha/2, n, p) + (-2:2)
    u[u > n] <- Inf
    l[l < 0] <- -Inf
    coverage <- outer(l, u, function(a,b) pbinom(b-1,n,p) - pbinom(a-1,n,p))
    if (max(coverage) < 1-alpha) i <- which(coverage==max(coverage)) else
      i <- which(coverage == min(coverage[coverage >= 1-alpha]))
    i <- i[1]
    # return the order statistics and the actual coverage.
    u <- rep(u, each=5)[i]
    l <- rep(l, 5)[i]
    return(list(Interval=c(l,u), Coverage=coverage[i]))
    }
  h = quantile.CI(reps(a),p,1-conf)
  if (h$Coverage < conf) cat('planned',conf,'actual',h$Coverage,'\n')
  sort(a@x)[h$Interval]
  }
## both Morgan & Henrion's qci() and Bill Huber’s quantile.CI function give confidence intervals for fractiles
## but they differ, and those differences become larger with small sample sizes
#savemany = MC$many
#MC$many = 1000
#conf = 0.90
##a = exponential(5)
##edf(a,new=TRUE,xlim=c(0,22))
#a = T(0,0,9)
#a = N(5,1)
#edf(a,new=TRUE)
##for (p in 2:3/10) {
##for (p in 1:9/10) {
#for (p in c(0.025,0.05,1:9/10,0.95,0.975)) {
#  q = qci(a,p, conf)
#  lines(c(q[[1]],q[[1]],q[[2]],q[[2]]), c(p,0,0,p))
#  lines(c(0,q[[2]]),c(p,p))
#  lines(rep(sort(a@x)[p*MC$many],2),c(0,p),col='red')
#  h = qCI(a,p,conf)   # sort(a@x)[quantile.CI(MC$many,p, 1-conf)$Interval]
#  lines(c(h[[1]],h[[1]],h[[2]],h[[2]]), c(p,0,0,p),col='green')
#  }
#MC$many = savemany
#
## the green intervals from Bill's quantile.CI() and the black intervals from Morgan & Henrion's qci() are a LITTLE DIFFERENT, but inconsistencies vary widely in different simulations with the same underlying distribution
## they seem to always agree about the left edges of intervals in the left tails (i.e., small p-values), and they commonly--but not not always--agree about the right edges in the right tails (this asymmetry itself seems strange)
## wrt the other bounds (right edges in the left tails, and left edges in the right tails), Bill's bounds are almost always shifted leftward in the right tail, and shifted rightward in the left tail compared to Morgan & Henrion
## I am not sure which we should prefer as it is possible they are both correct.  Their important differences may emerge when sample sizes are small.  Presumably we can compare them with COVERAGE SIMULATIONS.
#
#n <- 100      # Sample size
#q <- 0.90     # Percentile
##
## You only have to compute the order statistics once for any given (n,q).
##
#lu <- quantile.CI(n, q)$Interval
##
## Generate many random samples from a known distribution and compute 
## CIs from those samples.
##
#set.seed(17)
#n.sim <- 1e4
#index <- function(x, i) ifelse(i==Inf, Inf, ifelse(i==-Inf, -Inf, x[i]))
#sim <- replicate(n.sim, index(sort(rnorm(n)), lu))
##
## Compute the proportion of those intervals that cover the percentile.
##
#F.q <- qnorm(q)
#covers <- sim[1, ] <= F.q & F.q <= sim[2, ]
##
## Report the result.
##
#message("Simulation mean coverage was ", signif(mean(covers), 4), 
#        "; expected coverage is ", signif(quantile.CI(n,q)$Coverage, 4))


##########################################################################
# Wrap access methods
##########################################################################

#quiet <- setMethod('==', c('mc','mc'), function(e1, e2){ same(e1, e2) })

#if(!isGeneric("fivenum")) quiet <- setGeneric("fivenum", function (x, na.rm = FALSE) standardGeneric("fivenum"))
#quiet <- setMethod('fivenum', 'mc', function (x, na.rm = FALSE) fivenum.mc(x, na.rm = FALSE))
 
# there is no generic function iqr;  IQR returns the length of the iqr

if(!isGeneric("sd")) quiet <- setGeneric("sd", function(x, na.rm = FALSE) standardGeneric("sd"))
quiet <- setMethod('sd', 'mc', function(x, na.rm = FALSE) sd.mc(x))

if(!isGeneric("var")) quiet <- setGeneric("var", function(x, y = NULL, na.rm = FALSE, use) standardGeneric("var"))
quiet <- setMethod('var', 'mc', function(x, y, na.rm, use) var.mc(x))

if(!isGeneric("iqr")) quiet <- setGeneric("iqr", function(x, na.rm = FALSE) standardGeneric("iqr"))
quiet <- setMethod('iqr', 'mc', function(x) iqr.mc(x))

#this replaced the following hack, which appeared to work, but might collide with other packages 
#var <- function(x, ...) if (is.mc(x)) var.mc(x) else stats::var(x, ...)


##########################################################################
# Typing functions 
##########################################################################

is.mc <- function(x) inherits(x,'mc')

is.interval <- function(x) return(inherits(x,'interval') || (inherits(x,'pbox') && (x@d[[1]]==x@d[[steps(x)]]) &&  (x@u[[1]]==x@u[[steps(x)]])))

is.scalar <- function(x) {
  if (is.mc(x) && isTRUE(all.equal(left(x),right(x)))) return(TRUE)
  if (is.interval(x) && isTRUE(all.equal(left(x),right(x)))) return(TRUE)
# if (is.numeric(x) && identical(1,length(x))) return(TRUE)
  if (is.numeric(x) && isTRUE(1 == length(x))) return(TRUE)
  FALSE
  }

# whether it can be involved in a convolution
is.uncnum <- is.number <- function(x) is.numeric(x) || is.scalar(x) || is.interval(x) || is.mc(x) 

# whether it has epistemic or aleatory uncertainty
is.uncertain <- function(x) is.interval(x) || is.mc(x) 

is.vacuous <- function(s,...) {UseMethod("is.vacuous")};         is.vacuous.default    <- function(s) return(all(is.nan(s)) || all(is.na(s)))
is.vacuous.mc     <- function(s) all(s@u==-Inf) && all(s@d==Inf)
is.vacuous.interval <- function(s) all(s@lo==-Inf) && all(s@hi==Inf)
#is.vacuous.numeric  use default method

##########################################################################
# Unary transformations 
##########################################################################

abs.mc <- function(x) mc(abs(x@x))

exp.mc <- function(x) mc(exp(x@x))

log.mc <- function(x, expo=exp(1)) if (expo == exp(1)) mc(log(x@x)) else mc(log(x@x,expo))    # expo must be a scalar

log.mc <- function(x, expo=mc(exp(1))) mc(log(x@x,expo@x))    
  	
reciprocate.mc <- function(x) mc(1/x@x)

sqrt.mc <- function(x) mc(base::sqrt(x@x)) 
  
sin.mc <- function(x) mc(base::sin(x@x))   
cos.mc <- function(x) mc(base::cos(x@x))   
tan.mc <- function(x) mc(base::tan(x@x))   
asin.mc <- function(x) mc(base::asin(x@x))   
acos.mc <- function(x) mc(base::acos(x@x))   
atan.mc <- function(x) mc(base::atan(x@x))   
atan2.mc <- function(y,x) mc(base::atan2(y@x,x@x))   

`-.mc` <- negate.mc <- function(x) mc(-x@x) 

complement.mc <- function(x) mc(1-x@x) 

mult.mc <- function(m, x) if (m < 0) negate.mc(mult(abs(m),x)) else mc(m*x@x)

shift.mc <- function(ss,x) mc(ss+x@x) 

sign.mc <- function(x) mc(sign(x@x))

trunc.mc <- function(x) mc(base::trunc(x@x)) 
round.mc <- function(x,...) mc(base::round(x@x,...)) 
ceiling.mc <- function(x) mc(base::ceiling(x@x)) 

truncate.mc <- function(x,min,max) mc(base::pmin(max, base::pmax(min,x@x))) 

samedistribution <- function(x) mc(shuffle(x@x))

shuffle.mc <- function(x) mc(shuffle(x@x))


##########################################################################
# Mass reassignments functions 
##########################################################################

lowest <- function(x, p) {   # exclude all but the lowest p% (PERCENT, not fraction on [0,1]) of x, useful for Will Powley's rejection rescaling
  if (!is.mc(x)) stop('The argument to the lowest truncation must be a distribution')
  z <- sort(x@x)[1:(length(x@x)*(p/100))]
  mc(z[trunc(runif(MC$many,1,length(z)+1))])
  }

highest <- function(x, p) {   # exclude all but the highest p% (PERCENT, not fraction on [0,1]) of x, useful for Will Powley's rejection rescaling
  if (!is.mc(x)) stop('The argument to the highest truncation must be a distribution')  
  z = sort(x@x)[(length(x@x)*((100-p)/100)):length(x@x)]
  mc(z[trunc(runif(MC$many,1,length(z)+1))])
  }

below <- function(x, s) {   # exclude all values from x but those below or equal to s
  if (!is.mc(x)) stop('The argument to the below truncation must be a distribution')
  z = x@x[x@x <= s]
  mc(z[trunc(runif(MC$many,1,length(z)+1))])
  }

above <- function(x, s){   # exclude all values from x but those above or equal to s
  if (!is.mc(x)) stop('The argument to the above truncation must be a distribution')
  z = x@x[s <= x@x]
  mc(z[trunc(runif(MC$many,1,length(z)+1))])
  }

between <- function(x, m, M) below(above(x, m), M)   # exclude all values from x except those between m and M

rescale <- function(x, m, M)  m + (M-m) * (x - left(x))/ diff(range(x))  # linearly rescale x to the specified range [m,M]

##########################################################################
# K-fold binary convolutions, mixtures, selfsums, and extremizations
##########################################################################

# pmin, pmax are convolutions, corresponding to minI, maxI in Risk Calc

pmin.numeric = function(...,na.rm=FALSE){e=list(...); if (class(e[[2]])=='mc') return(pmin.mc(mc(e[[1]]),e[[2]],na.rm=na.rm)) else return(base::pmin(...,na.rm=na.rm))}
pmax.numeric = function(...,na.rm=FALSE){e=list(...); if (class(e[[2]])=='mc') return(pmax.mc(mc(e[[1]]),e[[2]],na.rm=na.rm)) else return(base::pmax(...,na.rm=na.rm))}

pmin.mc <- function (..., na.rm = FALSE) {  
  elts <- makemc(...)
  m <- elts[[1]]
  for (each in elts[-1]) m <- conv.mc(m, each, 'pmin')
  m
  }

pmax.mc <- function (..., na.rm = FALSE) {  
  elts <- makemc(...)
  m <- elts[[1]]
  for (each in elts[-1]) m <- conv.mc(m, each, 'pmax')
  m
  }

smax = function (..., na.rm = FALSE) {  
  elts <- makemc(...)
  x <- sort(elts[[1]]@x)
  for (each in elts[-1]) {
    x <- pmax(x,sort(each@x),na.rm=na.rm)
    }
  mc(shuffle(x))
  }
  
smin = function (..., na.rm = FALSE) {  
  elts <- makemc(...)
  x <- sort(elts[[1]]@x)
  for (each in elts[-1]) {
    x <- pmin(x,sort(each@x),na.rm=na.rm)
    }
  mc(shuffle(x))
  }
  
 mixture = function(x, w=rep(1,length(x)), ...) {   
  k <- length(x)
  if (k != length(w)) stop('Need same number of weights as arguments to mixture')
  w = w / sum(w)
  X = NULL
  for (i in 1:k)  X = c(X, mc(x[[i]]))
  x = X
  z = NULL
  for (i in 1:k)  z = c(z,random(x[[i]],w[[i]]*MC$many))
  z = shuffle(z)
  if (length(z)>MC$many) z = z[1:MC$many] else if (length(z)<MC$many) z = c(z, random(z, MC$many - length(z)))
  mc(z)
  }

extremize = extremise = function(f,k,max=TRUE) { 
	# assumes iid deviates
	# cf. pmax and pmin
	if (class(f)=='function') F = f else F = function() f
	if (max) {
		b = mc(-Inf)
		for (i in 1:k) b = pmax(b, shuffle(F())) 
		}
	else {
		b = mc(Inf) 
		for (i in 1:k) b = pmin(b, shuffle(F()))
		}
	return(b)
	}
#pl(0,10)
#gray(normal(5,1))
#cyan(extremise(normal(5,1),90))  # reuses the deviates
#salmon(extremise(function() normal(5,1), 900))  # uses fresh deviates for each max/min

selfsum = function(f,k) {   # search for other "selfsum" partial implementations
	b = mc(0)
	if (class(f)=='function') F = f else F = function() f
	for (i in 1:k) b = b + shuffle(F()) 
	return(b)
	}
#pl(0,10)
#gray(normal(5,1))
#cyan(selfsum(normal(5,1),20)/20)  # reuses the deviates
#salmon(selfsum(function() normal(5,1), 900)/900)  # uses fresh deviates for each addition


##########################################################################
# Binary convolution operations 
##########################################################################

conv.mc <- function(x,y,op='+') mc(do.call(op, list(x@x, y@x)))

perfectconv.mc <- function(x,y,op='+')  mc(shuffle(do.call(op, list(sort(x@x), sort(y@x)))))

oppositeconv.mc <- function(x,y,op='+')  mc(shuffle(do.call(op, list(sort(x@x), sort(y@x,decreasing=TRUE)))))


#####################
# correlation functions 
#####################

correlations <- function(S) {
  require('MASS'); 
  Sigma = matrix(S, nrow=sqrt(length(S))); 
  MC$r <- pnorm(mvrnorm(n=MC$many, mu=rep(0,nrow(Sigma)), Sigma=Sigma)); 
  return(invisible(NULL))}

R = function(i,corrs=MC$r) return(corrs[,i])

wranglenames = function(names) {
  names = unlist(strsplit(names, ','))
  names[[1]] = substring(names[[1]],3)
  nth = names[[length(names)]] 
  names[[length(names)]] = substring(nth, 1, nchar(nth)-1)
  names = trimws(names)
  }

corrs <- function(X,names=wranglenames(deparse(substitute(X)))) { 
  x = NULL
  k = length(X)
  for (i in 1:k) x = c(x,X[[i]]@x)
  C = cor(matrix(x,ncol=k))
  dimnames(C) <- list(names,names)
  C
  }
#corrs(c(A,B,C,D,E))

plotcorrs = function(X,names=wranglenames(deparse(substitute(X))),some=300,lwd=2,...) {
  x = NULL
  k = length(X)
  for (i in 1:k) x = c(x,X[[i]]@x)
  C = cor(matrix(x,ncol=k))
  par(mfcol=c(k,k))
  #w = shuffle(1:length(X@x))[1:some]
  w = shuffle(1:MC$many)[1:some]
  for (i in 1:k) for (j in 1:k) if (i==j) plot.mc(X[[i]],xlab=names[[i]],ylab='Cum.Prob.',lwd=lwd,...) else plot(X[[i]]@x[w],X[[j]]@x[w],xlab=names[[i]],ylab=names[[j]],main=paste(sep='','r=',signif(C[i,j],3)),...)
  invisible(NULL)
  }

# The caracel function will generate deviates with an exact specified sample (rather than population) correlation
# (This function will not work if the %*% matrix multiplication operation has been hijacked.  If the expression
# matrix(1:4,nrow=2) %*% matrix(3:6,nrow=2) doesn't give the 2-by-2 matrix (15   23; 22   34), you can use 
# rm("%*%") to restore its functionality, but doing this may break something else.)  The approach used by the 
# caracel function seems only to work with bivariate normal variables.  For instance, it does not work with 
# discrete random variable, because the trick used is rotation (try plotting x1 and x that result from letting 
# x2=rpois(2000,2); the resulting x2 is not Poisson).  This was caracal's answer on 
# http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
# See also http://www.r-bloggers.com/modelling-dependence-with-copulas-in-r/
# Use these lines to test the caracel function:
#    n = 2000
#    x1=rnorm(n,1,1)
#    x2=rpois(n, 2)
#    x2=rexp(n, 2)
#    x2=rnorm(n, 2, 0.5)
#    o = caracel(0.6, x1, x2)
#    X1 = o$x1
#    X2 = o$x2
#    cor(x1, X1); plot(sort(x1),1:n,type='s'); lines(sort(X1),1:n,type='s',col='red')
#    cor(x2, X2); lines(sort(x2),1:n,type='s',col='blue'); lines(sort(X2),1:n,type='s',col='tan')
#    cor(X1, X2)
caracel <- function(rho=0.6, x1=rnorm(2000,1,1), x2=rpois(2000, 2)) { # rho is the desired correlation = cos(angle)
  n = length(x1)                                               # n=20,000 almost crashed the computer!
  if (n!=length(x2)) stop('Mismatched vectors in caracel function')
  theta <- acos(rho)                                          # corresponding angle
  X     <- cbind(x1, x2)                                     # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  Id   <- diag(n)                                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))          # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)                                      # = Q Q'     # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                         # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                          # bind to matrix
  Y <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1   # needs the native R %*% matrix multiplication operator
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]             # final new vector
  x = ((x-mean(x))  / sd(x)) * sd(x2) + mean(x2)             # rescale to dispersion and location [added by Scott]
  list(x1=x1,x2=x,cor=cor(x1, x))                          # check correlation = rho
  }


#####################
# Logical functions 
#####################

# We don't support ==, &, &&, |, || as infix operators for p-boxes because 
# we cannot be comprehensive in doing so.  The first problem is that && 
# and || (which would be the natural operators to use) are sealed and 
# can't be changed.  The second problem is that & and I cannot be applied
# when the left argument is an atomic numeric for similar reasons.  Finally,
# & and | won't work on vectors of distributions, or rather, we cannot invoke 
# them that way, because an array of distributions does not have the class mc.
# Instead, we define the functions and, or, not, xor, nand, nor, cond, imply.

compare = function(a,b,tol=1e-10) {s = ifelse(abs(a-b)<tol, paste('Same:',a), paste('Different:' ,a,b)); cat(s,'\n'); invisible(s)}

are.logicals.numeric <- function(a,b) {
  # must be dimensionless too, but we cannot check that
  ((0 <= a) && (a <= 1) && (0 <= b) && (b <= 1))
  }

are.logicals.mc <- function(a,b) {
  a <- makemc(a)
  b <- makemc(b)
  # must be dimensionless too, but we cannot check that
  ((0 <= left.mc(a)) && (right.mc(a) <= 1) && (0 <= left.mc(b)) && (right.mc(b) <= 1))
  }

complement.numeric <- function(a) 1-a

	# if r=0, assumes independence
	# if r=+1, assumes perfect dependence
	# if r=-1, assumes opposite dependence
	# otherwise assumes Pearson correlation r, using the Lucas model, truncated as described in "Correlated Boolean Operators for Uncertainty Logic" by Miralles-Dolz et al.

and.numeric = function(p, q, r=0)  {
    # if (class(q)) return(and.mc(mc(p),q,r=r)) crashes for scalars
    Sandc = function(p , q , r ) { # scalar correlated conjunction
        denominator = sqrt(p*(1-p)*q*(1-q))
	if (abs(denominator)<1e-8) denominator = 1
        Lr = (pmax(p+q-1,0)-p*q)/denominator
        Ur = (pmin(p,q)-p*q)/denominator
        if (r < Lr) return(pmax(p+q-1,0))
        if (r > Ur) return(pmin(p,q))
        d = p * q + r*sqrt(p*q*(1-p)*(1-q))
        return(d)
        }
    lB = Sandc(p,q,r)
    uB = Sandc(p,q,r)
    if (1e-4<abs(lB-uB)) return(c(lB,uB)) else return(uB)
    }

or.numeric = function(p, q, r=0) 1 - and(1-p, 1-q, r)

	# if r=0, assumes independence
	# if r=+1, assumes perfect dependence
	# if r=-1, assumes opposite dependence
	# otherwise assumes Pearson correlation r, using the Lucas model, truncated as described in "Correlated Boolean Operators for Uncertainty Logic" by Miralles-Dolz et al.

and.mc = function(p, q, r=0)  {
    rightside.mc = leftside.mc = function(x) x
    Sandc = function(p , q , r ) { # scalar correlated conjunction
        denominator = sqrt(p*(1-p)*q*(1-q))
	if (abs(denominator)<1e-8) denominator = 1
        Lr = (pmax(p+q-1,0)-p*q)/denominator
        Ur = (pmin(p,q)-p*q)/denominator
        if (r < Lr) return(pmax(p+q-1,0))
        if (r > Ur) return(pmin(p,q))
        d = p * q + r*sqrt(p*q*(1-p)*(1-q))
        return(d)
        }
    # checkUncBool.([p,q]); checkCor(r)
    lB = Sandc(leftside.mc(p),leftside.mc(q), leftside.mc(r))
    uB = Sandc(rightside.mc(p),rightside.mc(q), rightside.mc(r))
    if (1e-4<abs(lB-uB)) return(c(lB,uB)) else return(uB)
    }

not = function(p) 1 - p
or = function(p, q, r=0) 1 - and(1-p, 1-q, r)
xor = function(p, q, r=0) p + q - 2 * and(p, q, r)
cond = function(p, q, r=0) and(p,q,r) / q
nand = function(p,q,r=0) complement(and(p, q, r))
nor = function(p,q,r=0) complement(or(p, q, r))
imply = function(p,q,r=0) or(complement(p), q, r)

# complement.mc (already defined) is logical negation

# and.mc <- function(a,b) conv.mc(mc(a),mc(b),'*'); # assumes independence!
#and.mc <- function(a,b,r=cor(a,b)) andc(mc(a),mc(b),r)

# or.mc <- function(a,b) complement.mc(conv.mc(complement.mc(mc(a)), complement.mc(mc(b)),'*')) # assumes independence!
#or.mc <- function(a,b,r=cor(a,b)) orc(mc(a),mc(b),r)

#cond.mc <- function(p,q,r=cor(a,b)) return(andc(p,q,r) / q)

# Dependence between arguments is automatic for arithmetic operations
# because the mc structures embody the dependence.  The same cannot 
# be said for logical operations because of what we call the "two-level-
# dependence issue").  By default, the mc functions use the correlation 
# between the operands to select the logical operation to use.  This may 
# or may not make sense.
#par(mfrow=c(4,2))
#MC$many=50000
#a = N(0.5, 0.08)
#b = U(.2,.4,r=a)
#cor(a,b)
#a
#b
#plot(a,b)
#aabi = andc(a,samedistribution(b),0)
#aab = andc(a,b)
#aabc = andc(a,b,cor(a,b))
#aabr = andc(a,b,+1)
##pl = function(...,xlab='',ylab='Cumulative probability') plot(NULL,ylim=c(0,1),xlim=range(...),xlab=xlab,ylab=ylab)
#pl(0,0.5,xlab='independent');   red(aabi)
#pl(0,0.5,xlab='default');           black(aab)
#pl(0,0.5,xlab='r=cor(a,b)');      green(aabc)
#pl(0,0.5,xlab='specified r=+1'); blue(aabr)
#pl(0,0.5,xlab='(blue is under green)'); black(aab); red(aabi); blue(aabr); green(aabc)


##########################################################################
# In-fix Operators #
##########################################################################

# addition 
quiet <- setMethod('+', c('mc','numeric'),  function(e1, e2){ shift.mc(e2, e1) })    # what happens if e2 is a vector??
quiet <- setMethod('+', c('numeric','mc'),  function(e1, e2){ shift.mc(e1, e2) })
quiet <- setMethod('+', c('mc','mc'),function(e1,e2){  return( conv.mc(e1,e2,'+') )})   

# multiplication
quiet <- setMethod('*', c('mc','numeric'), function(e1, e2){ mult.mc(e2, e1) })   # what happens if e2 is a vector??
quiet <- setMethod('*', c('numeric','mc'), function(e1, e2){ mult.mc(e1, e2) })
quiet <- setMethod('*', c('mc','mc'),function(e1,e2){ return( conv.mc(e1,e2,'*') )})   

# subtraction
quiet <- setMethod('-', c('mc','numeric'), function(e1, e2){mc(e1@x - e2) })   # what happens if e2 is a vector??
quiet <- setMethod('-', c('numeric','mc'), function(e1, e2){ mc(e1 + negate.mc(e2), opposite=e2) })
quiet <- setMethod('-', c('mc','mc'),function(e1,e2){ return( conv.mc(e1,negate.mc(e2),'+') )})   

# division
quiet <- setMethod('/', c('mc','numeric'), function(e1, e2){ mult.mc(1/e2, e1) })   # what happens if e2 is a vector??
quiet <- setMethod('/', c('numeric','mc'), function(e1, e2){ mult.mc(e1, reciprocate.mc(e2)) })
quiet <- setMethod('/', c('mc','mc'),function(e1,e2){ return( conv.mc(e1,reciprocate.mc(e2),'*') )})   

# raising to a power
quiet <- setMethod('^', c('mc','numeric'), function(e1, e2){ mc(e1@x ^ e2) })      # what happens if e2 is a vector??
quiet <- setMethod('^', c('numeric','mc'), function(e1, e2){ mc(e1 ^ e2@x)})
quiet <- setMethod('^', c('mc','mc'),function(e1,e2){ return( conv.mc(e1,e2,'^') )})   

MClt <- function(x,y) if (right(x) < left(y)) return(1) else if (right(y) < left(x)) return(0) else if (is.scalar(x) && is.scalar(y)) return(left(x)<left(y)) else return(prob.mc(conv.mc(mc(x),mc(y),'-'),0))
MCgt <- function(x,y) if (right(y) < left(x)) return(1) else if (right(x) < left(y)) return(0) else if (is.scalar(x) && is.scalar(y)) return(left(y)<left(x)) else return(prob.mc(conv.mc(mc(y),mc(x),'-'),0))
MClte <- function(x,y) if (right(x) <= left(y)) return(1) else if (right(y) < left(x)) return(0) else if (is.scalar(x) && is.scalar(y)) return(left(x)<=left(y)) else return(xprob.mc(conv.mc(mc(x),mc(y),'-'),0))
MCgte <- function(x,y) if (right(y) <= left(x)) return(1) else if (right(x) < left(y)) return(0) else if (is.scalar(x) && is.scalar(y)) return(left(y)<=left(x)) else return(xprob.mc(conv.mc(mc(y),mc(x),'-'),0))
MCeq <- function(x,y) {eq = mc(x)@x == mc(y)@x; return(sum(eq+0)/MC$many)}

inside.mc <- function(x,lo,hi=right(lo)) {hi=hi+0; lo=left(lo); length(x@x[(lo<x@x) & (x@x<hi)])/length(x@x)} # interestingly, the calculation "hi=hi+0" is needed to subvert lazy evaluation which would otherwise prevent inside.mc(N(5,1),c(0,5)) from getting the right answer

quiet <- setMethod('%in%', c('mc','numeric'),  function(x, table){ inside.mc(x, table) })   # N(5,1)  %in%  c(0,5)

# less-than comparison 
quiet <- setMethod('<', c('mc','numeric'), function(e1, e2){ MClt(e1,e2) })     # what happens if e2 is a vector??
quiet <- setMethod('<', c('numeric','mc'), function(e1, e2){ MClt(e1,e2) })
quiet <- setMethod('<', c('mc','mc'), function(e1, e2){ prob.mc(conv.mc(e1,negate.mc(e2),'+'),0) })

# greater-than comparison 
quiet <- setMethod('>', c('mc','numeric'),   function(e1, e2){ MCgt(e1,e2) })   # what happens if e2 is a vector??
quiet <- setMethod('>', c('numeric','mc'),  function(e1, e2){ MCgt(e1,e2) })
quiet <- setMethod('>', c('mc','mc'),  function(e1, e2){ xprob.mc(conv.mc(e2,negate.mc(e1),'+'),0) })

# less-than-or-equal-to comparison 
quiet <- setMethod('<=', c('mc','numeric'), function(e1, e2){ MClte(e1,e2) })     # what happens if e2 is a vector??
quiet <- setMethod('<=', c('numeric','mc'), function(e1, e2){ MClte(e1,e2) })
quiet <- setMethod('<=', c('mc','mc'), function(e1, e2){ xprob.mc(conv.mc(e1,negate.mc(e2),'+'),0) })

# greater-than-or-equal-to comparison 
quiet <- setMethod('>=', c('mc','numeric'),   function(e1, e2){ MCgte(e1,e2) })   # what happens if e2 is a vector??
quiet <- setMethod('>=', c('numeric','mc'),  function(e1, e2){ MCgte(e1,e2) })
quiet <- setMethod('>=', c('mc','mc'),  function(e1, e2){ prob.mc(conv.mc(e2,negate.mc(e1),'+'),0) })

# equality comparison (remember, for continuous distributions, equality comparions should be zero)
#a = N(5,1)
#a == 5           # 0
#a == a           #  1
#a == N(5,1)   # 0
#a == abs(a)   # 1
quiet <- setMethod('==', c('mc','numeric'), function(e1, e2){ MCeq(e1,e2) })     # what happens if e2 is a vector??
quiet <- setMethod('==', c('numeric','mc'), function(e1, e2){ MCeq(e1,e2) })
quiet <- setMethod('==', c('mc','mc'), function(e1, e2){ MCeq(e1,e2) })


##########################################################################
# Some inverse probability distributions not already implemented in R 
##########################################################################

qexponentialpower <- function(p, lambda, kappa) return(exp(-(log(lambda) - log(log(1-log(1-p)))) / kappa))

#qfrechet <- function(p, b, c) b*exp((-1/c)*log(log(1/p)))               # Castillo, page 207
qfrechet <- function(p, b, c, m=0) m+b*exp((-1/c)*log(log(1/p)))  

qgammaexponential <- function(p, shape, scale=1, rate=1/scale) rate * ((1-p)^(-1/shape) - 1)

qgev <- function(p,a,b,c) if (c==0) qgumbel(p,a,b) else a + b * (exp(-c*log(-log(p))) - 1) / c   # c must be a scalar

qgeneralizedpareto <- function(p,mu,sigma,scale) mu - sigma * (1-exp(-scale*log(1-p))) / scale

qgumbel <- function(p, a, b) a - b * log(log(1.0/p))                    # Evans et al., page 65

qlaplace <- function(p, a, b) ifelse(p<=0.5, a + b * log(2.0 * p),a - b * log(2.0 * (1.0 - p))) # Evans et al., page 92

qloglogistic <- function(p, lambda, kappa) return(exp(-(log(lambda) - log(-p/(p-1))/ kappa)))

qlogtriangular <- function(p, i_min, i_mode, i_max){
  a = log(i_min)
  b = log(i_mode)       # could this really be correct?
  c = log(i_max)
  exp(qtriangular(p, a, b, c))
  }

qloguniform <- function(p, one, two) {
  m = log(one)
  exp((log(two) - m) * p + m);
  }

qlomax <- function(p, lambda, kappa) return(((1-p)^(-1/kappa)-1)/lambda)

qreciprocal <- function(p, a=10.0) exp(p * log(a))

qpareto <- function(p, mode, c) mode * exp((-1.0/c) * log(1.0 - p))     # Evans et al., page 119

qpowerfunction <- function(p, b, c) b * exp((1.0/c) * log(p))           # Evans et al., page 128

qrayleigh <- function(p, b) sqrt(-2.0 * b * b * log(1.0 - p))           # Evans et al., page 134

sawinconradalpha01 <- function(mu) {
  if (abs(mu-0.5)<0.000001) return(0)
  f = function(alpha) 1/(1-1/exp(alpha)) - 1/alpha - mu
  uniroot(f,c(-500,500))$root
  }

qsawinconrad <- function(p, min, mu, max){
  alpha = sawinconradalpha01((mu-min)/(max-min))
  if (abs(alpha)<0.000001) return(min+(max-min)*p) else
  min+(max-min)*((log(1+p*(exp(alpha)-1)))/alpha)
  }
 
qshiftedloglogistic <- function(p, a,b,c) a + b * (exp(c*log(1/(1/p-1))) - 1) / c

qtrapezoidal <- function(p, a, b, c, d){
  if (abs(d-a)<1e-10) return(rep(a,length.out=length(p)))
  if (abs(c-b)<1e-10) return(qtriangular(p,a,b,d))
  h = 2 / (c+d-b-a)
  p1 = h * (b-a)/2
  p2 = p1 + h * (c-b)
  r <- ifelse(p <= p2, (p - p1) / h + b, d - sqrt(2 * (1-p) * (d-c) / h))
  r[p<=p1] <- a + sqrt(2 * p[p<=p1] * (b-a)/h)   
  r
  }

qtriangular <- function(p, min, mid, max){
  pm = (mid-min) / (max-min)
  q = ifelse(p<=pm, min + sqrt(p * (max - min) * (mid - min)), max - sqrt((1.0 - p) * (max - min) * (max - mid)))
  q[(max<mid) | (mid<min)] = NA
  q   
  }


##########################################################################
# Further distribution wrangling functions not already implemented in R 
##########################################################################

rtriangular <- function(many, min, mid, max)   if((max<mid) | (mid<min)) return(NA) else return(qtriangular(runif(many), min, mid, max))

# dtriangular <- function(x, min, mid, max) {ba = max - min; pmax(0,pmin(1,ifelse(x<mid, 2*(x-min)/ba*(mid - min), 2*(max-x)/(ba*(max - mid)))))}
dtriangular <- function(x, a, c, b) {ba=b-a; ca=ba*(c-a); bc=ba*(b-c); d=pmax(0,pmin(1,ifelse(x<c, 2*(x-a)/ca, 2*(b-x)/bc))); d[(b<c)|(c<a)]=NA; d}

#ptriangular <- function(x, a, c, b) {ba = b - a; ca = ba*(c - a); bc = ba*(b - c); ifelse(x<a,0,ifelse(b<x,1,ifelse(x<c, (x-a)^2/ca, 1-(b-x)^2/bc)))}
ptriangular <- function(x, a, c, b) {ba=b-a; ca=ba*(c-a); bc=ba*(b-c); p=ifelse(x<a,0,ifelse(b<x,1,ifelse(x<c, (x-a)^2/ca, 1-(b-x)^2/bc))); p[(b<c)|(c<a)]=NA; p}

#par(mfrow=c(5,2))
#x = -100:1000/100
#for (m in 0:9) {plot(x, ptriangular(x,0,m,7)); lines(x, dtriangular(x,0,m,7), col='red'); title(paste('0  ',m,'  7'))}
#plot(x, ptriangular(x,0,seq(0,10,length.out=length(x)),7))
#par(mfrow=c(4,2))
#x = -300:1000/100
#for (m in -2:5) {plot(x, ptriangular(x,-2,m,5),ylab=''); lines(x, dtriangular(x,-2,m,5), col='red'); title(paste('triangular(-2, ',m,', 5)'))}


##########################################################################
# Distribution constructors 
##########################################################################

#r_ <- function(r) if (!inherits(r,'mc')) return(r) else return(rank(r@x)/MC$many) 
#r_ <- function(r) if (!inherits(r,'mc')) return(r) else return(((1:MC$many)/(MC$many+1))[rank(r@x)])  # is this one slightly better?
#r_ = function(r) if (mode(r)=='list') {z1 = r[[2]]; rho = r[[1]]; return(pnorm(rho * (z1@x-mean(z1@x))/sd(z1@x) + sqrt(1-rho*rho)*rnorm(MC$many)))} else if (inherits(r,'mc')) return(rank(r@x)/MC$many) else return(r) 

#I think this versin of r_() should be slightly better, than the next one, but this one interacts badly with the dice() function
r_ = function(r) if (mode(r)=='list') {z1 = r[[2]]; rho = r[[1]]; return(pnorm(rho * (z1@x-mean(z1@x))/sd(z1@x) + sqrt(1-rho*rho)*rnorm(MC$many)))} else if (inherits(r,'mc')) return(((1:MC$many)/(MC$many+1))[rank(r@x)]) else return(r) 
#
r_ = function(r) if (mode(r)=='list') {z1 = r[[2]]; rho = r[[1]]; return(pnorm(rho * (z1@x-mean(z1@x))/sd(z1@x) + sqrt(1-rho*rho)*rnorm(MC$many)))} else if (inherits(r,'mc')) return(rank(r@x)/MC$many) else return(r) 
#
###############################
#rbyc(3,4)
#a = triangular(1,2,12)
#b = normal(10,2,r=c(0.75,a))
#c = normal(25,3, r=a)
#d = normal(25,3)
#plot(a);  title('a=T(1,2,12)')
#plot(b);  title('b=N(10,2,r=c(0.75,a))')
#plot(c);  title('c=N(25,3, r=a)')
#plot(d);  title('d=N(25,3)')
#plot(a,a); title(cor(a,a))
#plot(a,b); title(cor(a,b))
#plot(a,c); title(cor(a,c))
#plot(a,d); title(cor(a,d))
#e = -a
#f = normal(10,2,r=c(-0.75,a))
#g = normal(25,3, r=-a)
#h = samedistribution(a)
#plot(a,e); title(cor(a,e))
#plot(a,f); title(cor(a,f))
#plot(a,g); title(cor(a,g))
#plot(a,h); title(cor(a,h))
###############################
#rbyc(3,4)
#a = triangular(1,2,12)
#for (r in c(-0.99, -0.9, -0.75, -0.5, -0.3, -0.2, 0.2, 0.3, 0.5, 0.75, 0.9, 0.99)) {
#  b = normal(10,2,r=c(r,a))  
#  plot(a,b); title(r)
#  }
##############################

# fairdice() cannot use the dependence r argument, although dice() and coin() do
fairdice = function(sides=6) mc(sample(1:sides,MC$many,replace=TRUE))

dice = function(sides=6, p=rep(1/sides,sides)) mc(shuffle(rep(1:sides, ceiling(MC$many * (p/sum(p))))[1:MC$many])) 
dice = function(sides=6, p=rep(1/sides,sides),r) if (missing(r)) mc(shuffle(rep(1:sides, ceiling(MC$many * (p/sum(p))))[1:MC$many])) else mc(rep(1:sides, ceiling(MC$many * (p/sum(p)))) [MC$many*r_(r)])
#
# This example proves that it's possible to correlate positively through discrete distributions:
#####################
#MC$many = 1000
#a = N(5,1)
#R = (rep(1:sides, ceiling(MC$many * (p/sum(p)))))
#plot(a@x,R[MC$many*r_(a)])
#MC$many = 20000
#####################
#
# In fact, using  r_ = function(r) if (mode(r)=='list') {z1 = r[[2]]; rho = r[[1]]; return(pnorm(rho * (z1@x-mean(z1@x))/sd(z1@x) + sqrt(1-rho*rho)*rnorm(MC$many)))} else if (inherits(r,'mc')) return(rank(r@x)/MC$many) else return(r) 
# the dice() function with correlations seems to work beautifully:
#rbyc(2,4)
#p = 1:6
#a = N(5,1);           pl(a); red(a);   title('a=N(5,1)')
#b = dice(p=p);        pl(b); blue(b);  title('b=dice()')
#c = dice(p=p,r=a);    pl(c); green(c); title('c=dice(r=a)')
#d = dice(p=p,r=-a);   pl(d); purple(d);title('d=dice(r=-a)')
#plot(a,a); title(paste('cor(a,a) =',signif(cor(a,a),3)))
#plot(a,b); title(paste('cor(a,b) =',signif(cor(a,b),3)))
#plot(a,c); title(paste('cor(a,c) =',signif(cor(a,c),3)))
#plot(a,d); title(paste('cor(a,d) =',signif(cor(a,d),3)))
#
#
# The perfectconv.mc() and oppositeconv.mc() functions fail to preserve the natural 
# dependence between their inputs and the convolution result.  Thus it would be much 
# better if the dice() function were equipped to handle the dependence r parameter.
# HOWEVER, when we try to use the correlated dice() function to do opposite dependence
# IT LETS US DOWN.
#
#R = dice(p=1:6)
#B = dice(p=1:6)
######################### independence 
#cor(R,B)
#i = R * B
#pl(0,40); purple(i)
######################### perfect
#p = perfectconv.mc(R, B, '*')
#pl(0,40); orange(p)
#B = dice(p=1:6, r=R)
#cor(R,B)   # [1] 1
#pp = R*B
#brown(pp, lwd=2, lty='dashed')
########################## opposite
#o = oppositeconv.mc(R, B, '*')
#pl(0,40); gray(o)
#B = dice(p=1:6,r=-R)
#cor(R,B)
#oo = R*B
#black(oo,lwd=2,lty='dashed')
## oops!
#plot(R)     # MC (min=1, median=5, mean=4.33295, max=6) 
#plot(B)     # MC (min=2, median=4, mean=4.1909, max=6) 
## clearly B is corrupted, and not an opposite version of R
## or perhaps there cannot be an opposite version of R in 
## this discrete case 
#
## could another approach work with unfair dice??
#R = dice(p=1:6)
#B = dice(p=1:6)
#o = oppositeconv.mc(R, B, '*')
#pl(0,40); gray(o)
#B = 7 - R
#ooo = R*B
#black(ooo,lwd=2,lty='dashed')
#plot(R,B)
#cor(R,B)   # [1] -1
#  ??
#

#coin = function(p=0.5) mc((p < runif(MC$many))+0) # yields Heads (0) and Tails (1) where p is the probability of Heads
coin = function(p=0.5,r=runif(MC$many)) bernoulli(p,r) # yields Heads (1) and Tails (0) where p is the probability of Heads

bernoulli = function(p,r=runif(MC$many))  mc((r_(r) < p) + 0)
 
B <- beta  <- function(shape1, shape2,r=runif(MC$many)) mc(qbeta(r_(r), shape1, shape2))

beta1 <- function(mean,sd) return(beta(mean * (mean * (1 - mean) / (sd^2) - 1), (mean * (mean * (1 - mean) / (sd^2) - 1)) * (1/mean - 1)))

BB <- betabinomial <- function(size, v, w, r=runif(MC$many)) mc(qbinom(r_(r), size, beta(v,w)@x))  # pbox.r has a completely different implementation
 
Bin <- binomial <- function(size=NULL, prob=NULL, mean=NULL, std=NULL, r=runif(MC$many)){
  if (is.null(size) & !is.null(mean) & !is.null(std)) size <- mean/(1-std^2/mean)
  if (is.null(prob) & !is.null(mean) & !is.null(std)) prob <- 1-std^2/mean
  mc(qbinom(r_(r), size, prob))
  }

cauchy <- function(location, scale, r=runif(MC$many)) mc(qcauchy(r_(r), location, scale))

chi <- function(n, r=runif(MC$many)) mc(sqrt(chisquared(n)))

chisquared <- function(df, r=runif(MC$many)) mc(qchisq(r_(r),df))

delta <- dirac <- function(x, r=runif(MC$many)) mc(rep(x,MC$many))

discreteuniform <- function(max=NULL, mean=NULL, r=runif(MC$many)) {
  if (is.null(max) & !is.null(mean)) max <- 2 * mean
  mc(trunc(uniform(0, max+1)))
  }

exponential <- function(mean=1/rate, rate=1, r=runif(MC$many)) mc(qexp(r_(r),1/mean))

exponentialpower <- function(lambda, kappa, r=runif(MC$many)) mc(qexponentialpower(r_(r),lambda, kappa))
  
F <- f <- fishersnedecor <- function(df1, df2, r=runif(MC$many)) mc(qf(r_(r), df1, df2))

F1 = f1 = fishersnedecor1  = function(m,s) {Fw=2/(1-1/m); return(fishersnedecor((2*Fw^3  - 4*Fw^2)   /((Fw-2)^2 * (Fw-4) * s^2  - 2*Fw^2), Fw))}

#frechet <- function(b, c, r=runif(MC$many)) mc(qfrechet(r_(r), b, c))
#frechet <- function(scale, shape, r=runif(MC$many)) mc(qfrechet(r_(r), scale, shape))
frechet <- function(scale, shape, m=0, r=runif(MC$many)) mc(qfrechet(r_(r),scale,shape,m))

# N.B. the parameters for gamma are scale and shape, which is different from rgamma's parameters which are shape and rate (1/scale)
gamma <- function(scale=1/rate, shape, rate=1, r=runif(MC$many)) mc(qgamma(r_(r), shape=shape, scale=scale)) #mean: scale*shape, var: scale^2*shape

gamma2 <- function(shape, rate=1, scale=1/rate, r=runif(MC$many)) mc(qgamma(r_(r), shape=shape, scale=scale)) #mean: shape*scale, var: shape*scale^2 

gamma1 <- function(mu, sigma, r=runif(MC$many)) return(gamma(sigma^2/mu, (mu/sigma)^2, r=r_(r))) 

inversegamma <- function(shape, rate = 1, scale = 1/rate, r=runif(MC$many)) mc(rev(1/qgamma(r_(r), shape=shape, scale=scale))) # mean: 1/(scale*(shape-1)), var: 1/(rate^2*((shape-1)^2*(shape-2)))
  
# some older implementations
# gamma <- function(scale, shape, r=runif(MC$many)) mc(qgamma(r_(r), shape=shape, scale=scale)) 
# gamma1 <- function(mu, sigma) return(gamma(sigma^2/mu, (mu/sigma)^2)) 
# gamma2 <- function(shape, scale, r=runif(MC$many)) mc(qgamma(r_(r), shape=shape, scale=scale))  
# inversegamma <- function(a, b, r=runif(MC$many)) mc(1/qgamma(r_(r), shape=a, scale=1/b))  

gammaexponential <- function(shape, scale=1, rate=1/scale, r=runif(MC$many)) mc(qgammaexponential(r_(r),shape=shape,rate=rate))

geometric <- pascal <- function(prob=NULL, mean=NULL, r=runif(MC$many)){
  if (is.null(prob) & !is.null(mean)) prob <- 1/(1+ mean)
  mc(qgeom(r_(r),prob))
  }

GEV <- gev <- generalizedextremevalue <- fishertippett <- function(a=0,b=1,c=0, r=runif(MC$many)) if (c==0) return(gumbel(a,b)) else mc(qgev(r_(r),a,b,c)) # http://en.wikipedia.org/wiki/Generalized_extreme_value_distribution

generalizedpareto <- function(mu,sigma,scale, r=runif(MC$many)) mc(qgeneralizedpareto(r_(r),mu,sigma,scale))

gumbel <- extremevalue <- function(a=NULL, b=NULL, mean=NULL, std=NULL, var=NULL, r=runif(MC$many)){  
  if (missing(std) && !missing(var)) std <- sqrt(var)
  if (missing(a) && missing(b)) {a <- mean-std*0.577215665*sqrt(6)/base::pi; b <- std*sqrt(6)/base::pi}
  mc(qgumbel(r_(r), a, b))
  }

laplace <- function(a, b, r=runif(MC$many)) mc(qlaplace(r_(r), a, b))  

logistic <- function(location, scale, r=runif(MC$many)) mc(qlogis(r_(r),location, scale))

loglogistic <- fisk <- function(lambda, kappa, r=runif(MC$many)) mc(qloglogistic(r_(r),lambda, kappa))

lognormal0 <- function(meanlog, stdlog, r=runif(MC$many)) mc(qlnorm(r_(r),meanlog,stdlog)) 

L <- lognormal <- function(mean=NULL, std=NULL, meanlog=NULL, stdlog=NULL, median=NULL, cv=NULL, r=runif(MC$many)){
  if (is.null(meanlog) & !is.null(median)) meanlog = log(median)
  if (is.null(stdlog) & !is.null(cv)) stdlog = sqrt(log(cv^2 + 1))
  # lognormal(a, b) ~ lognormal2(log(a^2/sqrt(a^2+b^2)),sqrt(log((a^2+b^2)/a^2)))
  if (is.null(meanlog) & (!is.null(mean)) & (!is.null(std))) meanlog = log(mean^2/sqrt(mean^2+std^2))
  if (is.null(stdlog) & !is.null(mean) & !is.null(std)) stdlog = sqrt(log((mean^2+std^2)/mean^2))
  if (!is.null(meanlog) & !is.null(stdlog)) lognormal0(meanlog,stdlog,r) else stop('Not enough information to specify the lognormal distribution')
  }

logtriangular <- function(min=NULL, mode=NULL, max=NULL, minlog=NULL, midlog=NULL, maxlog=NULL, r=runif(MC$many)){
  if (is.null(min) & !is.null(minlog)) min = exp(minlog)  ###place in arglist
  if (is.null(max) & !is.null(maxlog)) max = exp(maxlog)  ###place in arglist
  mc(qlogtriangular(r_(r), min, mode, max))
  }

# loguniform1 needs loguniform_f, LUgrid and loguniform_solve 
loguniform_f = function(a,m,v) a*m*exp(2*(v/(m^2)+1)) + exp(2*a/m)*(a*m - 2*((m^2) + v))
LUgrid = function(aa, w) left(aa)+(right(aa)-left(aa))*w/100.0   # aa is an interval
loguniform_solve = function(m, v) {
  aa = c(m - sqrt(4*v), m)   # interval
  a = m
  ss = loguniform_f(a,m,v)
  for (j in 1:4) {
    for (i in 0:100) {
      a = LUgrid( aa, i)
      s = abs(loguniform_f(a,m,v))
      if (s<ss) {ss = s; si = i }
    }
    a = LUgrid(aa, si)
    aa = c(LUgrid(aa, si-1), LUgrid(aa, si+1))  # interval
  }
  return(a)
}

loguniform <- function(min=NULL, max=NULL, minlog=NULL, maxlog=NULL, mean=NULL, std=NULL, r=runif(MC$many)){
  if (is.null(min) & !is.null(minlog)) min <- exp(minlog)
  if (is.null(max) & !is.null(maxlog)) max <- exp(maxlog)
  if (is.null(max) &  !is.null(mean) & !is.null(std) & !is.null(min)) max = 2*(mean^2 +std^2)/mean - min
  if (is.null(min) & is.null(max) &  !is.null(mean) & !is.null(std)) {
    min = loguniform_solve(mean,std^2)
    max = 2*(mean^2 +std^2)/mean - min
    }
  mc(qloguniform(r_(r), min, max))
  }

loguniform1 = function(m, s) loguniform(mean=m, std=s)

lomax <- function(lambda, kappa, r=runif(MC$many)) mc(qlomax(r_(r),lambda, kappa))   # support = (0,infinity)      (scale) lambda > 0, (shape) kappa > 0

NB <- negativebinomial <- function(size, prob, mu, lower.tail = TRUE, log.p = FALSE, r=runif(MC$many)) mc(qnbinom(r_(r), size, prob, mu, lower.tail = lower.tail, log.p = log.p ))
#rr = c(1,2,3,4,5,10,20,40)
#kk = 0:25
#for (r in rr) {
#  p = 1/(1+r/10)  
#  plot(kk, dnbinom(kk,r,p))
#  for (k in kk) lines(rep(k,2),c(0,dnbinom(k,r,p)),col='cornflowerblue',lwd=4)
#  title(paste('r=',r))
#  Sys.sleep(2)
#  }

normal0 <- function(normmean, normstd, r=runif(MC$many)) mc(qnorm(r_(r),normmean,normstd))

N <- normal <- gaussian <- function(mean=NULL, std=NULL, median=NULL, mode=NULL, cv=NULL, IQR=NULL, var=NULL, px=NULL, r=runif(MC$many)){
  if (is.null(mean) & !is.null(median)) mean = median
  if (is.null(mean) & !is.null(mode)) mean = mode
  if (is.null(std) & !is.null(var)) std = sqrt(var)
  if (is.null(std) & !is.null(cv) & !is.null(mean)) std = mean * cv
  if (is.null(std) & !is.null(IQR)) std = IQR / 1.34896    
  if (!is.null(mean)) px = c(px,c(0.5,mean))
  if ((is.null(mean) || is.null(std)) & !is.null(px)) { # px specifies any two percentiles, e.g., px = c(0.5, 0, 0.84134, 1) specifies the standard normal
    q = qnorm(px[[3]]); 
    if (q == 0) {
        if (is.null(mean)) mean = px[[4]]
	std = (px[[2]] - mean) / qnorm(px[[1]])
        } 
    else {
        Q = qnorm(px[[1]]) / q
	if (is.null(mean)) mean = (px[[2]] - Q*px[[4]]) / (1 - Q)
        std = (px[[4]] - mean) / q
        }
    }
  if (!is.null(mean) & !is.null(std)) normal0(mean,std,r) else stop('Not enough information to specify the normal distribution')
  }

sawinconrad <- function(min, mu, max, r=runif(MC$many)){
  a <- left(min);   b <- right(max)
  c <- left(mu);    d <- right(mu)
  if (c<a) c <- a   # implicit constraints
  if (b<d) d <- b
  mc(qsawinconrad(r_(r), min, mu, max))
  }

SN <- skewnormal <- function(location, scale, skew, r=runif(MC$many)) if (require('sn')) mc(qsn(r_(r),location,scale,skew)) else stop('Need to install the sn package to use skewnormal')

pareto <- function(mode, c, r=runif(MC$many))  mc(qpareto(r_(r), mode, c))

poisson <- function(lambda, r=runif(MC$many)) mc(qpois(r_(r),lambda))

poissonbinomial <- function(p, r=runif(MC$many)) {  # p is an array of (possibly different) probability values
  # the distribution of the sum of independent Bernoulli trials, with possibly different probabilities of success p
  # support = 0:(length(p)), mean = sum(p), variance = sum(p(1-p))
  # poissonbinomial(p[1],p[2],p[3],...,p[n]) ~ binomial(n, p*) iff p[i]=p*, i=1,...n
  # http://en.wikipedia.org/wiki/Poisson_binomial_distribution, 8 May 2012
  # Wang, Y. H. (1993). "On the number of successes in independent trials". Statistica Sinica 3 (2): 295-312.
  n = length(p)
  m = sum(p)
  v = sum(p*(1-p))
  Prk = rep(0,n+1)
  Prk[0+1] = prod(1-p)
  T <- function(i) sum((1/(1/p-1))^i)    #sum((p[j]/(1-p[j]))^i)
  Ti = rep(0,n)
  for (i in 1:n) Ti[i] = T(i)
  for (k in 1:n) {
    i = 1:k
    Prk[k +1] = (1/k) * sum((-1)^(i-1) * Prk[k-i +1] * Ti[i])
    }
  C = cumsum(Prk)
  u = rep(0,MC$many)
  i = r_(r)
  for (k in 1:n) u[C[k] <= i] = k
  mc(u)
  }

powerfunction <- function(b, c, r=runif(MC$many)) mc(qpowerfunction(r_(r), b, c))

rayleigh <- function(b, r=runif(MC$many)) mc(qrayleigh(r_(r), b)) 

reciprocal <- function(a=10.0, r=runif(MC$many)) mc(qreciprocal(r_(r), a))

shiftedloglogistic <- function(a=0,b=1,c=0, r=runif(MC$many))  {
  if (c==0) return(Slogistic(a,b))
  csc <- function(x) 1/sin(x)
  m <- a + b * (base::pi*c*csc(base::pi*c)) / c
  v <- b^2  * (2*base::pi*c*csc(2*base::pi*c)-(base::pi*c*csc(base::pi*c))^2) / (c^2)
  mc(qshiftedloglogistic(r_(r),a,b,c))
  }

student <- function(df, r=runif(MC$many)) mc(qt(r_(r),df))

trapezoidal <- function(min, lmode, rmode, max, r=runif(MC$many)) mc(qtrapezoidal(r_(r), min, lmode, rmode, max))
  
T <- triangular <- function(Min, Mode, Max, r=runif(MC$many)) mc(qtriangular(r_(r), Min, min(Mode,Max), max(Mode,Max)))
  
U <- uniform <- rectangular <- # function(min, max, r=runif(MC$many)) mc(qunif(r_(r),min,max))
function(min=NULL,max=NULL,mean=NULL,sd=NULL,median=NULL,cv=NULL,r=runif(MC$many)) {
  if (!is.null(min) & !is.null(max)) return(mc(qunif(r_(r),min,max)))
  if (is.null(mean) & !is.null(median)) mean = median
  if (is.null(sd) & !is.null(cv) & !is.null(mean)) sd = mean * cv
  if (!is.null(mean) & !is.null(sd)) return(mc(qunif(r_(r),mean-sqrt(3)*sd,max=mean+sqrt(3)*sd))) 
  else stop('Not enough information to specify the uniform distribution')
  }
    
# argument order disagrees with R's qweibull function, but agrees with Wikipedia (in November 2013)
weibull <- function(scale, shape, r=runif(MC$many)) mc(qweibull(r_(r), shape, scale))   


##########################################################################
# Custom distribution constructors 
##########################################################################

histogram = function(x) { #sample values are in an array x 
  n <- length(x)
  i <- shuffle(rep(1:n, ceiling(MC$many/n)))[1:MC$many]
  # k <- rep(0,n); for (j in 1:n) k[j] <- length(i[i==j]); if (length(i)!=MC$many) cat(n,' ',length(i),' ',k,'\n')
  mc(x[i])
  }
  
quantiles <- function(v,p,r=runif(MC$many))  {
  if (length(v) != length(p)) stop('Inconsistent array lengths for quantiles')
  if ((min(p) < 0) || (1 < max(p))) stop('Improper probability for quantiles') # ensure 0 <= p <= 1
  if (!identical(range(p),c(0,1))) stop('Probabilities must start at zero and go to one for quantiles')
  if (any(diff(p)<0)) stop('Probabilities must increase for quantiles') # ensure montone probabilities
  if (any(diff(v)<0)) stop('Quantiles values must increase') # ensure montone quantiles
  r = r_(r)
  x = rep(Inf,MC$many)
  for (i in 1:(length(p)-1)) x = ifelse((p[[i]]<=r) & (r<p[[i+1]]), v[[i]]+(r-p[[i]])*(v[[i+1]]-v[[i]])/(p[[i+1]]-p[[i]]),  x)
  mc(x)
  }
#quantiles(v = c(3, 6,   10,  16,  25,  40), p = c(0, 0.2, 0.4, 0.5, 0.7, 1))
#quantiles(v = c(3,   6,   10,  16,  25,  40), p = c(0.1, 0.2, 0.4, 0.5, 0.7, 0.9))


##########################################################################
# PERT and MaxEnt distribution constructors 
##########################################################################

betapert = function(min, max, mode) {
  mu = (min + max + 4*mode)/6
  if (abs(mode-mu)<1e-8) alpha1 = alpha2 = 3 else {
  alpha1 = (mu - min)*(2*mode - min - max)/((mode - mu)*(max - min))
  alpha2 = alpha1*(max - mu)/(mu - min) }
  min + (max - min) * beta(alpha1, alpha2)
  }
  
MEminmax <- function(min, max) uniform(min,max)

MEminmaxmean <- function(min, max, mean) sawinconrad(min,mean,max) #http://mathoverflow.net/questions/116667/whats-the-maximum-entropy-probability-distribution-given-bounds-a-b-and-mean, http://www.math.uconn.edu/~kconrad/blurbs/analysis/entropypost.pdf for discussion of this solution.

MEmeansd = function(mean, sd) normal(mean, sd)

MEminmean<- function(min,mean) min+exponential(mean-min)

#MEmeangeomean <- function(mean, geomean)

MEdiscretemean <- function(x,mu,steps=10,iterations=50) { # e.g., MEdiscretemean(1:10,2.3)
  fixc = function(x,r) return(1/sum(r^x))
  r = br = 1
  c = bc = fixc(x,r)
  d = bd = (mu - sum((c*r^x)*x))^2
  for (j in 1:steps) {
    step = 1/j
    for (i in 1:iterations) {
      r = abs(br + (runif(1) - 0.5) * step)
      c = fixc(x,r)
      d = (mu - sum((c*r^x)*x))^2
      if (d < bd) {
        br = r
        bc = c
        bd = d
        }
      }
    }
  w = bc*br^x
  w = w / sum(w) # needed?
  z <- NULL
  k = length(x)
  for (i in 1:k)  z <- c(z,rep(x[[i]],w[[i]]*MC$many))
  if (length(z)>MC$many) z = z[1:MC$many] else if (length(z)<MC$many) z = c(z, z[shuffle(z)[1:(MC$many - length(z))]])
  mc(shuffle(z))
  }

MEquantiles <- quantiles

MEdiscreteminmax <- function(min,max) return(pmin(trunc(uniform(min,max+1)),max))

MEmeanvar <- function(mean, var) return(MEmeansd(mean,sqrt(var)))

MEminmaxmeansd <- function(min, max, mean, sd) return(beta1((mean- min) / (max - min),  sd/(max - min) ) * (max - min) + min)

MEmmms <- MEminmaxmeansd

MEminmaxmeanvar <- function(min, max, mean, var) return(MEminmaxmeansd(min,max,mean,sqrt(var)))


##########################################################################
# Method-of-Moment distribution constructors (match moments of the data x)
##########################################################################

MMbernoulli <- function(x) bernoulli(mean(x)) # assumes x is zeros and ones
MMbeta <- function(x) beta1(mean(x), sd(x))

MMbetabinomial  = function(n,x) {
  # n must be provided; it's not estimated from data
  # https://en.wikipedia.org/wiki/Beta-binomial_distribution#Example:
  # MMbetabinomial(n=12,rep(0:12,c(3,24,104,286,670,1033,1343,1112,829,478,181,45,7))) 
  m1 = mean(x)
  m2 = mean(x^2)
  d = n*(m2/m1 - m1 - 1) + m1
  betabinomial(n, (n*m1 - m2) / d, (n-m1)*(n - m2/m1) / d)
  }
  
MMbinomial <- function(x) {a = mean(x); b= sd(x); binomial(round(a/(1-b^2/a)), 1-b^2/a)} 
MMchisquared <- function(x) chisquared(round(mean(x)))
MMexponential <- function(x) exponential(mean(x))
MMF <- function(x) {w = 2/(1-1/mean(x)); F(round((2*w^3 - 4*w^2) / ((w-2)^2 * (w-4) * sd(x)^2 - 2*w^2)), round(w))}
MMgamma <- function(x) {a = mean(x); b= sd(x); gamma(b^2/a, (a/b)^2)}  #gamma1(a, b) ~ gamma(b²/a, (a/b)²)
MMgeometric <- MMpascal <- function(x) geometric(1/(1+mean(x)))
MMgumbel <- MMextremevalue <- function(x) gumbel(mean(x) - 0.57721* sd(x) * sqrt(6)/ pi, sd(x) * sqrt(6)/ pi)
MMlognormal <- function(x) lognormal(mean(x), sd(x))
MMlaplace <- MMdoubleexponential <- function(x) laplace(mean(x), sd(x)/sqrt(2))
MMlogistic <-  function(x) logistic(mean(x), sd(x) * sqrt(3)/pi)
MMloguniform <- function(x) loguniform1(mean(x), sd(x))   
MMnormal <- MMgaussian <- function(x) normal(mean(x), sd(x))
MMpareto <- function(x) {a = mean(x); b= sd(x); pareto(a/(1+1/sqrt(1+a^2/b^2)), 1+sqrt(1+a^2/b^2))}
MMpoisson <- function(x) poisson(mean(x))
MMpowerfunction <- function(x) {a = mean(x); b= sd(x); powerfunction(a/(1-1/sqrt(1+(a/b)^2)), sqrt(1+(a/b)^2)-1)}
MMt <- MMstudent <- function(x) if (1<sd(x)) student(2/(1-1/sd(x)^2)) else stop('Improper standard deviation for student distribution')
MMuniform <- MMrectangular <- function(x) {a = mean(x); b= sd(x); uniform(a-sqrt(3)*b, a+sqrt(3)*b)}

MMuniform1 = function(w) {mu1 = mean(w)
mu2 = mean(w^2)
m = sqrt(3*(mu2-mu1^2))
uniform(mu1 - m, mu1 + m)}
# the example in https://en.wikipedia.org/wiki/Method_of_moments_(statistics)#Examples suggests
# the formula to use should be MMuniform1, but it doesn't yield results that match the data in moments:
#MC$many = 1000000
#x = runif(5)
#a = MMuniform(x)
#b = MMuniform1(x)
#plot(c(0,1),c(0,1),type='l',col='gray')
#black(a)
#blue(b)
#mean(x); mean(a); mean(b)
#sd(x); sd(a); sd(b)

MMtriangular = function(x,iters=100,dives=10) {
  # iterative search for triangular distribution parameters using method of matching moments (you solve the thing analytically! too messy without help)
  # testing code indicated with ##
  ##some = 10
  ##A = runif(1,0,10)
  ##B = A + runif(1,0,10)
  ##C = runif(1,A,B)
  ##x = qtriangular(runif(some), A,C,B)
  skewness = function(x) {m = mean(x); sum((x-m)^3)/((length(x)-1)*sd(x)^3)}
  M = mean(x)
  V = var(x)
  S = skewness(x)
  a = aa = min(x)
  b = bb = max(x)
  c = cc = 3*M-a-b 
  many = iters
  s1 = sd(x)
  for (k in 1:dives) {
    s1 = s2 = s3 = s1/2
    a = rnorm(many,aa,s1)
    b = rnorm(many,bb,s2)
    c = rnorm(many,cc,s3)
    m = (a+b+c)/3
    k = (a^2+b^2+c^2-a*b-a*c-b*c)
    v = k/18
    s = (sqrt(2)*(a+b-2*c)*(2*a-b-c)*(a-2*b+c)) / (5 * k ^ (3/2))
    d = (M-m)^2 + (V-v)^2 + (S-s)^2
    i = which.min(d)
    aa = a[[i]]
    bb = b[[i]]
    cc = c[[i]]
    }
  ##gray(triangular(A,B,C), new = TRUE)
  ##blue(x)
  ##green(triangular(aa,bb,cc))
  ##A;aa; B;bb; C;cc  # the order is min, max, mode
  return(triangular(a,c,b)) 
  }

triangular.antweiler = function(x, r=runif(MC$many)) triangular(Min=min(x), Mode=3*mean(x)-max(x)-min(x), Max=max(x), r=r) # https://wernerantweiler.ca/blog.php?item=2019-06-05

## triangular.antweiler respects the observed range of the data, that is, it finds a triangular distribution that does not go beyond this range
#tt = rtriangular(10,1,2,3)
#pl(0,4)
#edf(tt)
#blue(MMtriangular(tt))
#red(triangular.antweiler(tt))



#
#
#
#
#
#
#
#
#
#
#find_weibull = function(m,s,lambda_explore=NULL, k_explore=NULL) {
#  if (missing(lambda_explore)) lambda_explore = runif(50,m*0.5,m*1.5)
#  if (missing(k_explore)) k_explore = runif(50,0,400)
#
#  # initial estimates
#  lambda = m  
#  k = 400/(s*1.559483)
#
#
#  nl = nk = wm = wv = NULL
#  for (lambda in lambda_explore) for (k in k_explore) {
#    nl = c(nl, lambda)
#    nk = c(nk, k)
#    g = Gamma(1+1/k)
#    wm = c(wm, lambda * g)
#    wv =c(wv, (lambda^2) * (Gamma(1+2/k) - g^2))
#    #w = weibull(nl,nk)
#    }
#    
#  plot(wm,nl)
#  plot(wv,nl)
#  rl = lm(nl ~ wm + wv)
#
#  plot(wm,nk)
#  plot(wv,nk)
#  rk = lm(nk ~ wm + wv)
#  # list(rl=rl, rk=rk)
#  v = s^2
#  cat('Fitting Weibull to have mean & sd:', m, ' ', s,'\n')
#  lambda = predict(rl, newdata=data.frame(wm = m, wv = v))
#  k = predict(rk, newdata=data.frame(wm = m, wv = v))
#  cat('Fitted Weibull to have scale & shape:', lambda, ' ', k,'\n')
#  return(weibull(scale=lambda, shape=k))
#  }
#  
#w = find_weibull( 199.4277, 1.274248)
#mean(w); sd(w)
#
#
#
#
#
#lambda = 200
#k = 200
#
#g = Gamma(1+1/k)
#lambda * g                   # [1] 199.4277
#(lambda^2) * (Gamma(1+2/k) - g^2)   # [1] 1.623708
#sqrt((lambda^2) * (Gamma(1+2/k) - g^2))  # [1] 1.274248
#
#w = weibull(lambda,k)
#mean(w)   # [1] 199.4429
#var(w)    # [1] 1.550305
#sd(w)    # [1] 1.245112
#
#
#
#lambda = 100
#ss =NULL
#for (k in 1:400) {
#  g = Gamma(1+1/k)
#  s = sqrt((lambda^2) * (Gamma(1+2/k) - g^2))  # [1] 1.274248
#  ss = c(ss,s)
#  cat(lambda * g, '\t',                   # [1] 199.4277
#  #(lambda^2) * (Gamma(1+2/k) - g^2)   # [1] 1.623708
#  s
#  , '\n')
#  }
#  
#plot(1:400,ss)
#
#plot(1:400,1/ss)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#k = 200:2000/1000
#G2 = Gamma(1+2/k)
#plot(G2,k)
#R = lm(k ~ poly(G2, 10))
#K = predict(R, newdata=data.frame(G2=G2))
#lines(G2,K,col=2)
#cor(k,K)
#plot(k,K)
#
#IG2 = function(g2) predict(R, newdata=data.frame(G2=g2))
#
#MMweibull = function(x) {
#  m = mean(x)
#  s = sd(x)
#  v = s^2
#  # initial estimates
#  lambda = m  
#  k = 400/(s*1.559483)
#cat('0 :', lambda, ' ', k,'\n')   
#  for (i in 1:20) {
#    g = Gamma(1+1/k)
#    lambda = m / g  
#    #Gamma(1+2/k)  = v / (lambda^2) + g^2
#    k = IG2(v / (lambda^2) + g^2)
#  cat(i,':', lambda, ' ', k,'\n')   
#    }
#  cat('Fitting Weibull to have mean & sd:', m, ' ', s,'\n')
#  cat('Fitted Weibull has scale & shape:', lambda, ' ', k,'\n')
#  return(weibull(scale=lambda, shape=k))
#  }
#
#
#W = weibull(200,10)
#pl(0,W)
#blue(W)
#data = rweibull(20, scale=200,shape=10)
#edf(data)
#w = MMweibull(data)
#red(w)
#edf(data)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#



##########################################################################
# Maximum likelihood estimation constructors
##########################################################################

MLbernoulli <- function(x) bernoulli(mean(x))
MLnormal <- MLgaussian <- MMnormal
MLexponential <- MMexponential
MLpoisson <- MMpoisson
MLgeometric <- MLpascal <- MMgeometric
MLuniform <- MLrectangular <- function(x) uniform(min(x), max(x))
MLpareto <- function(x) pareto(min(x), length(x)/sum(log(x)-log(min(x))))
MLlaplace <- MLdoubleexponential <- function(x) laplace(median(x), sum(abs(x-median(x))/length(x)))
MLlognormal_ = function(x) {n = length(x); mu = sum(log(x))/n; lognormal(meanlog=mu, stdlog=sum((log(x)-mu)^2)/n)}  # this function gives clearly poor results
MLlognormal = function(x) return(exp(MLnormal(log(x)))) # just uses transformation, which seems unlikely to be true, but fitdistrplus package uses it too
MLweibull <- function(x, shapeinterval=c(0.001,500)) { 
  f = function(k) sum(x^k * log(x)) / sum(x^k) - sum(log(x)) / length(x) - 1/k  
  k = uniroot(f, shapeinterval)$root
  el = exp(log(sum(x^k)/length(x))/k)
  weibull(scale=el, shape=k)
  }
MLgamma = function(data) {
  xbar = mean(data)
  shape=(xbar/sd(data))^2  # initial estimate of shape from MoM
  logxbar = log(xbar)
  meanlog = mean(log(data))
  f <- function(x) log(x) - digamma(x) - logxbar + meanlog
  shape = uniroot(f,shape*c(0.5,1.5))$root
  rate = shape/xbar
  return(gamma(shape=shape, rate=rate))
  }


##########################################################################
# Additional compound and conjugate distribution constructors for use in Bayesian inference
##########################################################################

# When there is no data, the posterior is the prior distribution.  
# As the sample size increases, the posterior tends to the data 
# manifested as the likelihood.  Typically, as the prior hyper-
# parmeters (a,b) increase, the posterior grows away from the 
# likelihood and towards the prior.  Generally, the default values 
# for the prior hyperparameters yield the uninformative case. 

# compound distributions already defined
# BB, betabinomial, gammaexponential, NB, negativebinomial, poissonbinomial

BCbernoulli <- function(x,a=0.5,b=0.5,only=TRUE) {
  # data x is an array of zeros and ones (failures and successes)
  # a and b are the parameters of the prior beta
  # see also the km( ) and KN( ) functions
  s = sum(x)
  n = length(x)
  lk = beta(s+1, n-s+1)
  pp = betabinomial(1, s+a, n-s+b)
  if (only) return(pp) else return(list(pr = beta(a, b), po = beta(s+a, n-s+b), pp = pp, lk = lk ))}
##############################################################################    
#par(mfrow=c(3,4))
#options(digits=3)
#x = c(0,0,0,1,1)
#  BCbernoulli(x,only=FALSE) #Jeffreys is the default
#  b = BCbernoulli(x,0,0)  # Haldane prior
#  b; title(paste('Haldane', mean(b)))
#  abline(h=1-sum(x)/length(x),col='green')
#  b = BCbernoulli(x)  # Jeffreys prior
#  b; title(paste('Jeffreys prior', mean(b)))
#  b = BCbernoulli(x,1,1)  # Bayes-Laplace prior
#  b; title(paste('Bayes-Laplace', mean(b)))
#  b = BCbernoulli(x,2,2)  # Walley prior
#  b; title(paste('Walley', mean(b)))
#x = c(1,1,rep(0,8))
#  BCbernoulli(x,only=FALSE) #Jeffreys is the default
#  b = BCbernoulli(x,0,0)  # Haldane prior
#  b; title(paste('Haldane', mean(b)))
#  abline(h=1-sum(x)/length(x),col='green')
#  b = BCbernoulli(x)  # Jeffreys prior
#  b; title(paste('Jeffreys prior', mean(b)))
#  b = BCbernoulli(x,1,1)  # Bayes-Laplace prior
#  b; title(paste('Bayes-Laplace', mean(b)))
#  b = BCbernoulli(x,2,2)  # Walley prior
#  b; title(paste('Walley', mean(b)))
#b = BCbernoulli(x,0,0,only=FALSE); Haldane = b$pr; pl(0,1, xlab='Haldane'); green(Haldane); title('prior')
#b = BCbernoulli(x,only=FALSE);      Jeffreys = b$pr;  pl(0,1, xlab='Jeffreys'); blue(Jeffreys); title('prior')
#b = BCbernoulli(x,1,1,only=FALSE); Laplace = b$pr;  pl(0,1, xlab='Laplace'); black(Laplace); title('prior')
#b = BCbernoulli(x,2,2,only=FALSE); Walley = b$pr;   pl(0,1, xlab='Walley'); gray(Walley); title('prior')

BCbinomial <- function(N, k,n,a=0.5,b=0.5,only=TRUE) {
  # data k is the count of successes, and n is a corresponding number of trials
  # both k and n may be arrays, but they must have the same length
  # a and b are the parameters of the prior beta
  # N is the number of trials to use for the posterior predictive distribution
  # see also the km( ) and KN( ) functions
  s = sum(k)
  sn = sum(n)
  lk = beta(s+1, sn-s+1)
  pp = betabinomial(N, s+a, sn-s+b)
  if (only) return(pp) else return(list(
    pr = beta(a, b),
    po = beta(s+a, sn-s+b),   
    pp = pp,
    lk = lk
    ))}
##############################################################################    
## https://stats.stackexchange.com/questions/512148/beta-binomial-vs-updating-a-prior-beta-distribution
## the BCbernoulli function is specified in terms of 1's (successes) and 0's (failures)
## the BCbinomial function is specified in terms of successes and number of TRIALS
## the beta and betabinomial distributions are specified in terms of numbers of successes and number of FAILURES 
#par(mfrow=c(1,1))
#b = BCbinomial(1,7,7+18,5,12,only=FALSE)
#bb = BCbernoulli(c(rep(1,7),rep(0,18)),5,12,only=FALSE)
#b
#bb
#pl(0,1)
#cyan(beta(5,12));  blue(b$pr);  cyan(bb$pr,lty='dotted')# should all be the same
#black(beta(7+5, 18+12));  gray(b$po);  black(bb$po,lty='dotted')# should all be the same
#khaki(betabinomial(1,12,30));  green(b$pp);  khaki(bb$pp,lty='dotted')# should all be the same
#pl(0,1)
#bc = BCbinomial(1,c(),c(),2,2,only=FALSE)
#blue(bc$pr)
#gray(bc$po)
#green(bc$pp)
##############################################################################    
#par(mfrow=c(3,5))
#options(digits=3)
#N = 10
#k = c(0,2)
#n = c(3,2)  # that is, the Bernoulli trials were {{0,0,0},{1,1}}
#2/5
#sum(k)/sum(n)
#  BCbinomial(N,k,n,only=FALSE) #Jeffreys is the default
#  b1 = BCbinomial(N,k,n,0,0)  # Haldane prior
#  pl(0,N); green(b1); title(paste('Haldane', mean(b1)))
#  b2 = BCbinomial(N,k,n)  # Jeffreys prior
#  pl(0,N); blue(b2); title(paste('Jeffreys prior', mean(b2)))
#  b3 = BCbinomial(N,k,n,1,1)  # Bayes-Laplace prior
#  pl(0,N); black(b3); title(paste('Bayes-Laplace', mean(b3)))
#  b4 = BCbinomial(N,k,n,2,2)  # Walley prior
#  pl(0,N); gray(b4); title(paste('Walley', mean(b4)))
#  pl(0,N); green(b1); blue(b2); black(b3); gray(b4)
## x = c(1,1,rep(0,8)) = {1,1,0,0,0,0,0,0,0,0}
#k = c(2,0,0,0)
#n = c(3,3,2,2)
#2/10
#sum(k)/sum(n)
#  BCbinomial(N,k,n,only=FALSE) #Jeffreys is the default
#  b1 = BCbinomial(N,k,n,0,0)  # Haldane prior
#  pl(0,N); green(b1); title(paste('Haldane', mean(b1)))
#  b2 = BCbinomial(N,k,n)  # Jeffreys prior
#  pl(0,N); blue(b2); title(paste('Jeffreys prior', mean(b2)))
#  b3 = BCbinomialN(k,n,1,1)  # Bayes-Laplace prior
#  pl(0,N); black(b3); title(paste('Bayes-Laplace', mean(b3)))
#  b4 = BCbinomial(N,k,n,2,2)  # Walley prior
#  pl(0,N); gray(b4); title(paste('Walley', mean(b4)))
#  pl(0,N); green(b1); blue(b2); black(b3); gray(b4)
#b1 = BCbinomial(N,k,n,0,0,only=FALSE); Haldane = b1$pr; pl(0,1, xlab='Haldane'); green(Haldane); title('prior')
#b2 = BCbinomial(N,k,n,only=FALSE);      Jeffreys = b2$pr;  pl(0,1, xlab='Jeffreys'); blue(Jeffreys); title('prior')
#b3 = BCbinomial(N,k,n,1,1,only=FALSE); Laplace = b3$pr;  pl(0,1, xlab='Laplace'); black(Laplace); title('prior')
#b4 = BCbinomial(N,k,n,2,2,only=FALSE); Walley = b4$pr;   pl(0,1, xlab='Walley'); gray(Walley); title('prior')
#pl(0,1); green(Haldane); blue(Jeffreys); black(Laplace); gray(Walley); title('prior')

BCpoisson = function(x,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the default hyperparameters seem to be the uninformative case
  s = sum(x)
  n = length(x)
  lk = gamma(shape=s+1, rate=n)
  pr = gamma(shape=a, rate=b)   #prr = rgamma(many,shape=a,rate=b)
  po = gamma(shape = a + s, rate = b + n)  
  #por = rgamma(r_(r),shape = a + s, rate = b + n)
  #ppr = mc(rpois(MC$many,por))   
  pp = negativebinomial(a + s, 1-1/(1+b + n))  
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}
##############################################################################    
#doit = function(x,a,b) {
#  pl(0,18)
#  points(x,rep(-0.017,length(x)),col='red')
#  bc = BCpoisson(x,a,b,only=FALSE)
#  blue(bc$pr)
#  gray(bc$po)
#  green(bc$pp)
#  edf(bc$ppr)
#  }
#par(mfcol=c(4,2))
#a = 2; b = 4
#doit(rpois(3000,5),a,b)   # should be the same as poisson(5)
#doit(rpois(10,5),a,b)   
#doit(c(3,4,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior
#a = 1; b = 1/50
#doit(rpois(3000,5),a,b)   # should be the same as poisson(5)
#doit(rpois(10,5),a,b)   
#doit(c(3,4,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior

BCgeometric = function(x,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the default values for the hyperparameters a and b are the uninformative case, although they will make the function crash if x=NULL
  s = sum(x)
  n = length(x)
  lk = beta(n+1, s+1)
  pr = beta(a,b)  
  po = beta(a + n, b + s)  
  por = rbeta(r_(r), a + n, b + s) 
  pp = mc(rgeom(MC$many,por))   # do we know an analytical formula?  
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}
##############################################################################    
#doit = function(x,a,b) {
#  pl(0,18)
#  points(x,rep(-0.017,length(x)),col='red')
#  bc = BCgeometric(x,a,b,only=FALSE)
#  blue(bc$pr)
#  gray(bc$po)
#  green(bc$pp)
#  }
#par(mfcol=c(4,2))
#a = 1/2; b = 1/2
#doit(rgeom(3000,0.5),a,b)   # should be the same as geometric(0.5)
#doit(rgeom(10,0.5),a,b)   
#doit(c(3,14,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior
#a = 1; b = 10
#doit(rgeom(3000,0.5),a,b)   # should be the same as geometric(0.5)
#doit(rgeom(10,0.5),a,b)   
#doit(c(3,14,1),a,b)
#doit(c(),a,b)                    # the posterior should equal the prior

BCuniform.knownmin = BCuniform = function(x,A,a=A,b=A+1,r=runif(MC$many),only=TRUE) {  
  #  x_i ~ uniform(A,theta), that is, from a uniform distribution whose minimum is A and whose maximum needs to be established
  lk = pareto(max(x), length(x)+1)
  pr = pareto(a, b)
  po = pareto(max(a,max(x)), b+length(x)) # Masatoshi says max(x) is m
  por = qpareto(r_(r),max(a,max(x)), b+length(x)) # Masatoshi says max(x) is m
  pp = mc(runif(MC$many,A,por))
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}
#x = runif(25,5,13)
#bc = BCuniform(x,A=5,a=1,b=1,only=FALSE) 
#bc

BCnegativebinomial = function(x,R,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the default hyperparameters are the uninformative case
  s = sum(x)
  n = length(x)
  lk = beta(R*n+1, s+1)
  pr = beta(a, b)
  po = beta(a+R*n, b+s)
  por = rbeta(r_(r),a+R*n, b+s)
  pp = mc(rnbinom(MC$many,R,por))   # do we know an analytical formula?  
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}

BCnormal.knownsigma = function(x,sigma,m0=mean(x),s0=sd(x),r=runif(MC$many),only=TRUE) {
  # the default hyperparameters seem to be uninformative (using s0=sd(x)/sqrt(length(x) makes the posterior differ from the likelihood more strongly)
  # increasing s0 makes the prior more uninformative, unlike the typical behaviour of (a,b) hyperparameters in other functions
  s = sum(x)
  n = length(x)
  pr = normal(m0,s0)  
  lk = normal(s/n,sigma/sqrt(n))
  mprime = (m0/s0^2+s/sigma^2)/(1/s0^2+n/sigma^2)
  sprime = 1/sqrt(1/s0^2+n/sigma^2)
  if (abs(s0<1e-20)) {mprime = m0; sprime = 0}
  po = normal(mprime, sprime)
  por = rnorm(r_(r),mprime, sprime)
  pp = normal(mprime, sqrt(sprime^2 + sigma^2))    # pp = mc(rnorm(MC$many,por,sigma)) 
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}

BCnormal.knownmu = function(x,mu,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the hyperparameter defaults are the uninformative case (so the pr won't be defined as it would be improper)
  n = length(x)
  s = sum((x-mu)^2)
  pr = sqrt(inversegamma(shape=a,rate=b)) # we parameterize N with sd, not var
  po = sqrt(inversegamma(shape=a + n/2, rate=b + s/2))
  por = sqrt(1/rgamma(MC$many,shape=a + n/2, rate=b + s/2))
  pp = mc(rnorm(MC$many,mu,por))     
  #ppa = mean(x)+sd(x)*sqrt(1+1/n)*student(n-1)   
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp))}
#  
#x = rnorm(1000,15,2)
#bc = BCnormal.knownmu(x,15,5,5)
#bc
#edf(x)
#  
#x = rnorm(10,15,2)
#bc = BCnormal.knownmu(x,15,5,5)
#bc
#edf(x)
#  
#x = rnorm(1000,15,2)
#bc = BCnormal.knownsigma(x,2,10,2)
#bc
#edf(x)
#  
#x = rnorm(10,15,2)
#bc = BCnormal.knownsigma(x,2,10,2)
#bc
#edf(x)

rnormgamma <- function(n, mu, lambda, alpha, beta) {
  # normal-gamma deviates: (1) Sample tau from a gamma distribution with parameters alpha and beta, (2) Sample x from a normal distribution with mean mu and variance 1/(lambda * tau)
  # E(x) = mu;  E(tau) = alpha/beta
  if (length(n) > 1) n <- length(n)
  tau <- rgamma(n, alpha, beta)
  x <- rnorm(n, mu, sqrt(1/(lambda*tau)))    # tau and x are NOT independent, as tau is used to compute x
  #cat(mean(x), mu,'  ', mean(tau), alpha/beta,'\n')
  data.frame(x = x, tau = tau)                              
}

seepr = function(m,l,a,b) {pr = rnormgamma(MC$many,m, l, a, b); prx = mc(pr$x); prt = mc(pr$t); pl(min(left(prx),left(prt)), max(right(prx),right(prt))); blue(prx); cyan(prt); title(paste(m,l,a,b)) }

BCnormal = function(x, mu0=mean(x), lambda0=1, alpha0=1, beta0=1, only=TRUE) {
  # x_i ~ N(mu, 1/tau)
  # The selection of the prior for a normal involves chosing values for 4 parameters.  The first argument 
  # mu0 is your guess about the mean of the normal data, and second is related to the dispersion of this 
  # estimate about the mean.  You can use the seepr( ) function to visualize the priors you select for 
  # BCnormal.  In practice, setting mu0 to mean(x) and the other values to one seems to often give 
  # reasonable results, so they have been specified as defaults, but I'm sure this violates some crucial
  # Bayesian stricture I'm not aware.
  pr = rnormgamma(MC$many,mu0, lambda0, alpha0, beta0)  # hmm...I guess pr and po and pp should all be correlated, right?  And they should have an intended correlation with ab outside caller via r_(r)
  n = length(x)
  xbar = mean(x)
  s = sd(x)
  po = rnormgamma(MC$many, (lambda0*mu0 + n*xbar)/(lambda0+n), lambda0+n, alpha0+n/2, beta0+(n*s+(lambda0*n*(xbar-mu0)^2)/(lambda0 + n))/2 )
  pp = mc(rnorm(MC$many,po$x,1/po$tau))
  pr = list(x=mc(pr$x), tau=mc(pr$tau))
  po = list(x=mc(po$x), tau=mc(po$tau))
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp))}
#x = rnorm(20,15000,2)
#m = mean(x)
#l = 1
#a = 1
#b = 1
#seepr(m,l,a,b)
#bc = BCnormal(x,m,l,a,b)
#bc
#edf(x)
#mean(bc)
#sd(bc)

BCexponential = function(x,a=0,b=0,r=runif(MC$many),only=TRUE) {
  # the prior and posterior estimate the MEAN of the exponential
  # which is the reciprocal of its RATE parameter used by rexp()
  # the default hyperparameters are the uninformative case
  sm = 'exponential(theta)'
  s = sum(x)
  n = length(x)
  lk = gamma(shape=n+1,rate=s)         # reciprocated in the returned list
  pr = gamma(shape=a, rate=b)           # reciprocated in the returned list
  po = gamma(shape=a+n, rate=b+s)  # reciprocated in the returned list
  por = rgamma(r_(r), shape=a+n,  rate=b+s)
  pp = mc(rexp(MC$many, por))
  if (only) return(pp) else return(list(pr = 1/pr, po = 1/po, pp = pp, lk = 1/lk))}

BCpareto.knownmin = BCpareto = function(x,xm,a=0,b=0,r=runif(MC$many),only=TRUE) { 
  # default hyperparameters seem to be the uninformative case
  s = sum(log(x/xm))
  n = length(x)
  sm = 'pareto(xm,theta)'
  lk = gamma(shape=n, rate=s)  # this is just a guess; prolly wrong as po not always between pr and this lk
  #pr = gamma(a, b)
  #po = gamma(1/(1/a+s), b+n) # po = gamma(1/(a+s), b+length(n))
  #por = rgamma(r_(r), 1/(1/a+s), b+n)
  #pp = mc(qpareto(r_(r),xm,por))
  pr = gamma(shape=a, rate=b)
  po = gamma(shape=a+n, rate=b+s) 
  por = rgamma(r_(r), shape=a+n, rate=b+s)
  pp = mc(qpareto(r_(r),xm,por))
  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp, lk = lk))}

#BCweibull <- function(x,a,b,c,d0,v=0) {
#  # x_i ~ weibull(shape=beta, scale=theta)
#  # https://www.johndcook.com/CompendiumOfConjugatePriors.pdf, page 33f
#  
#  n = length(x)
#  #sxb = sum(x^beta)
#  px = prod(x)
#  #lk = (beta/theta)^n * prod(x) ^ (beta-1) * exp(- sum(x^beta) / theta)
#  
#  #pr = function(beta,theta,a=a,b=b,c=c,d=d0,v=0) {K = 1;  D = function(beta,v,d) sum(d^beta[1:(v+1)]); ifelse((0<beta) & (0<theta), beta^(a-1)*exp(-beta*b)*theta^(-c)*exp(-D(beta,v,c(d0,x))/theta) /K,0)}
#  pr = function(beta,theta,a=a,b=b,c=c,d=d0,v=0) {K = 1;  D = function(beta,v,d) sum(d^beta[1:(v+1)]); beta^(a-1)*exp(-beta*b)*theta^(-c)*exp(-D(beta,v,c(d0,x))/theta) /K}
#  #need to evaluate normalization factor K
#  
#  po = function(beta,theta) pr(beta,theta, a=a+n, b=b-log(px), c=c+n, d0=d0, v=n) 
#
#  # the pr and po are bivariate, but the pp is just univariate, so it seems like we should be able to smash [not marginalize, right?] all the pr and po values into a single array with which to create the pp
#  # or should pobeta and potheta just be the marginalizations?
#
#  pp = mc(qweibull(r_(r),shape=pobeta,scale=potheta))
#  if (only) return(pp) else return(list(pr = pr, po = po, pp = pp))}

#x = qweibull(runif(30), shape=2, scale=3)
#x = c(0.9, 1.52, 1.10)
#a = 20.0
#b = 2.0
#c = 6.0
#d0 = 2.5
#v=0
#  n = length(x)
#  px = prod(x)
#pr = function(beta,theta,a=a,b=b,c=c,d=d0,v=0) {K = 1;  D = function(beta,v,d) sum(d^beta[1:(v+1)]); beta^(a-1)*exp(-beta*b)*theta^(-c)*exp(-D(beta,v,c(d0,x))/theta) /K} #need to evaluate normalization factor K
#po = function(beta,theta) pr(beta,theta, a=a+n, b=b-log(px), c=c+n, d0=d0, v=n) 
  
# old naming scheme  
#bernoulli.uniform = function(x) {s = sum(x); beta(1+s,1+length(x)-s)}
#bernoulli.beta = function(x,a,b) {s = sum(x); beta(a+s, b+length(x)-s)
#binomial.beta = function(x,N,a,b) {s = sum(x); beta(a+s, b+sum(N)-s)
#negativebinomial.beta = function(x,r,a,b)  beta(a+r*length(x), b+sum(x))
#poisson.gamma = function(x,k,theta) gamma(k + sum(x), 1/(length(x)+1/theta))
#poisson.gamma = function(x,a,b) gamma(a + sum(x), b + length(x))
#hypergeometric.betabinomial = function(x,N,a,b) betabinomial(a+sum(x), b+sum(N) - sum(x))
#geometric.beta = function(x,a,b) beta(a+length(x), b + sum(x))
#normal.normal = function(s, mu, sigma) {v=1/(sigma^2) + length(x)/s; normal((mu/(sigma^2) + sum(x)/(s^2)/v, 1/sqrt(v))) }
#exponential.gamma = function(x,a,b) gamma(a+length(x),b+sum(x))

km <- function(k,m) {
  # Bayesian posterior using the Jeffreys prior for the binomial rate 
  # (that is, the probability of success) given k successes and m failures 
  # randomly observed in k + m independent Bernoulli trials 
  if ((k < 0)  || (m < 0)) stop('Improper arguments to function km')
  return(beta(k+0.5,m+0.5))
  }

KN = function(k,n) {
  # Bayesian posterior using the Jeffreys prior for the binomial rate 
  # (that is, the probability of success) given only k successes out 
  # of n randomly observed independent Bernoulli trials 
  if ((k < 0)  || (n < k)) stop('Improper arguments to function KN')
  return(beta(k+0.5,n-k+0.5))
  }

# bc double-dot functions below don't have sampling models, just the likelihood functions

bc..beta <- function(x,c,d,a,b) {
  list(
    lk = beta(c,d),
    pr = beta(a,b),
    po = beta(a+c-1, b+d-1),
    pp = ''
    )}

bc..gamma <- function(x,d,e,b,c) {
  list(
    lk = gamma(d,e),
    pr = gamma(b,c),
    po = gamma(1/(1/b+1/d),c+e-1),
    pp = ''
    )}

bc..uniform <- function(x,c,d,a,b) {
  list(
    po = 0,
    lk = uniform(c,d),
    pr = uniform(a,b),
    po = uniform(max(a,c), min(b,d)),
    pp = ''
    )}

bc..normal <- function(x,c,d,a,b) {     # caution: b,d are repeated
  list(
    lk = normal(c,d),
    pr = normal(a,b),
    po = normal((a/(b^2)+c/(d^2))/(1/(b^2)+1/(d^2)),sqrt(1/(1/(b^2)+1/(d^2)))),
    pp = ''
    )}

bc..pareto <- function(x,c,d,a,b) {
  list(
    lk = pareto(c,d),
    pr = pareto(a,b),
    po = pareto(max(a,c), b+d+1),
    pp = ''
    )}

bc..powerfunction <- function(x,d,e,b,c) {
  list(
    lk = powerfunction(d,e),
    pr = powerfunction(b,c),
    po = powerfunction(min(b,d), c+e-1),
    pp = ''
    )}

bc..exponential <- function(x,b,a) {
  list(
    lk = exponential(b),
    pr = exponential(a),
    po = exponential(1/(1/a+1/b)),
    pp = ''
    )}

bc..laplace <- function(x,a,d,b) {
  list(
    lk = laplace(a,d),
    pr = laplace(a,b),
    po = laplace(a,1/(1/b+1/d)),      # lk and pr means must be the same,
    pp = ''
    )}

bc..gamma <- function(x) {
  list(
    lk = gamma(b,c),
    pr = exponential(a),
    po = gamma(1/(1/a+1/b),c),
    pp = ''
    )}

bc..exponential <- function(x) {
  list(
    lk = exponential(a),
    pr = gamma(b,c),
    po = gamma(1/(1/a+1/b),c),
    pp = ''
    )}

bc..paretobeta <- function(x) {
  list(
    lk = pareto(a,c),
    pr = beta(v,w),
    po = beta(v-c-1,w),
    pp = ''
    )}

bc..dirichlet <- function(x,a,b,r=runif(MC$many),only=TRUE) {
  list(
    lk = 0,
    pr = dirichlet(s, tj),
    po = dirichlet(n+s,(xj+stj)/(n+s)),
    pp = '' 
    )}


##########################################################################
# Odds language
##########################################################################

probability = function(wins=NULL, losses=NULL, of=wins/losses, oa=1/of) return(1/(1 + oa))

oddpair = function(p, lots=101) { # not good for extreme probabilities
  L = 1:lots;               A = rep(L,lots);   B = rep(L,rep(lots,lots));   P = A / (A+B)
  o = order(P);           P = (P[o]);          A = (A[o]);                      B = (B[o])
  d = !duplicated(P);   P = P[d];            A = A[d];                         B = B[d]
  i = which.min(abs(P-p))
  return(c(A[[i]],B[[i]]))
  }  

odds = function(p, how=1, lots=101, v=':') {# ratio of the prob it happens to the prob it doesn't: p / (1-p)
  of = 1 / (1/p - 1)          # odds in favor   
  if (how<0) of = 1/of     # odds against
  if (how == 1)  return(of)      # ? fractional, UK traditional  NOT SURE ABOUT THIS LINE
  if (how == -1) return(of)      # fractional, UK traditional      
  if (how == 2)  return(of+1)  # ? decimal, European  NOT SURE ABOUT THIS LINE
  if (how == -2) return(of+1)  # decimal, European
  if (how == 3)  if (of>= 1) return(100*of) else return(-100/of)  # ? American, Moneyline  NOT SURE ABOUT THIS LINE
  if (how == -3) if (of>= 1) return(100*of) else return(-100/of)  # American, Moneyline
  if (how == 4)  return(oddpair(p,lots))
  if (how == -4) return(rev(oddpair(p,lots)))
  if (how == 5)  return(paste(oddpair(p,lots),collapse=v))
  if (how == -5) return(paste(rev(oddpair(p,lots)),collapse=v))
  stop("Don't know how to compute those odds")
  }
#
#odds_infavor = odds_for = function(p,how,lots=100) odds(p,abs(how),lots)
#odds_against = function(p,how,lots=100) odds(p,-abs(how),lots)
#
#probability(1,1)
#probability(1,3)
#probability(1,4)
#odds(0.5)
#odds(0.2)
#odds(0.25)
#odds(0.99)
#rbyc(3,1)
#ppp = NULL
#pp = seq(0,1,length.out=100)
#for (p in pp) {o = odds(p,4,lots=100); ppp=c(ppp, probability(wins=o[[1]], losses=o[[2]]))}
#plot(pp,ppp,type='l'); lines(pp,pp,col=2)
#ppp = NULL
#pp = seq(0,1,length.out=100)
#for (p in pp) {o = odds(p,4,lots=10); ppp=c(ppp, probability(wins=o[[1]], losses=o[[2]]))}
#plot(pp,ppp,type='l'); lines(pp,pp,col=2)
#ppp = NULL
#pp = 0.25 + seq(-0.001,0.002,length.out=100)
#for (p in pp) {o = odds(p,4,lots=100); ppp=c(ppp, probability(wins=o[[1]], losses=o[[2]]))}
#plot(pp,ppp,type='l'); lines(pp,pp,col=2)
## https://www.aceodds.com/bet-calculator/odds-converter.html
#pp = c(.99,.833,.818,.8,.778,.769,.75,.733,.714,.692,.667,.652,.636,.619,.60,.579,.556,.545,.524,.5,.488,.476,.465,.455,.444,.421,.417,.4,.385,.381,.364,.357,.348,.333,.312,.308,.294,.286,.278,.267,.25,.238,.222,.20,.182,.167,.154,.143,.133,.125,.118,.111,.1,.091,.083,.077,.071,.067,.062,.059,.053,.048,.038,.029,.02,.015,.01,.001)
#tab = '\t'
#for (p in pp) cat(odds(p, -1), tab, odds(p,-5), tab, odds(p,-2), tab, odds(p,-3), tab, p,'\n')


##########################################################################
# Diagnostic Bayes' rule
##########################################################################

br = function(p,s,c) 1 / (1 + ((1/p - 1) * (1 - c)) / s)  


##########################################################################
# Decision making
##########################################################################

EU = function(m, p, rowsacts=TRUE) { # maximize expected utility for payoff matrix m and probabilities p;  choices are rows if rowsacts is TRUE
  if (!rowsacts) return(EU(base::t(m)))
  n = nrow(m)
  mm = rep(0,n)
  for (i in 1:n) mm[[i]] = sum(p * m[i,])  # how to vectorize?          expect = function(x) p*x; mm = apply(m,1,expect); cat(mm,nl)  
  M = max(mm)
  w = M == mm    
  r = dimnames(m)[[1]][w]
  if (is.null(r)) return(w) else return(r)
  }

maximin = function(m, rowsacts=TRUE) {
  if (!rowsacts) return(maximin(base::t(m)))
  n = nrow(m)
  mm = apply(m,1,min)  # mm = rep(0,n); for (i in 1:n) mm[[i]] = min(m[i,])  
  M = max(mm)
  w = M == mm    
  r = dimnames(m)[[1]][w]
  if (is.null(r)) return(w) else return(r)
  }

maximax = function(m, rowsacts=TRUE) {
  if (!rowsacts) return(maximax(base::t(m)))
  n = nrow(m)
  mm = apply(m,1,max)  
  M = max(mm)
  w = M == mm    
  r = dimnames(m)[[1]][w]
  if (is.null(r)) return(w) else return(r)
  }

hurwicz = function(m, h=0.5, rowsacts=TRUE) {
  if (!rowsacts) return(hurwicz(base::t(m)))
  n = nrow(m)
  mm = rep(0,n)
  for (i in 1:n) mm[[i]] = h*min(m[i,]) + (1-h)*max(m[i,])  # how to vectorize?
  M = max(mm)
  w = M == mm    
  r = dimnames(m)[[1]][w]
  if (is.null(r)) return(w) else return(r)
  }

bayes_laplace = chris_rock = function(m, rowsacts=TRUE) {
  if (!rowsacts) return(bayes_laplace(base::t(m)))
  n = nrow(m)
  mm = apply(m,1,sum) 
  M = max(mm)
  w = M == mm    
  r = dimnames(m)[[1]][w]
  if (is.null(r)) return(w) else return(r)
  }
  
minimaxregret = function(m, rowsacts=TRUE) {
  if (!rowsacts) return(minimaxregret(base::t(m)))
  n = nrow(m)
  mm = apply(m,2,max) # minuends
  regret = m
  for (i in 1:n) regret[i,] = mm - regret[i,] 
  mm = apply(regret,1,max)  
  M = min(mm)
  w = M == mm    
  r = dimnames(m)[[1]][w]
  if (is.null(r)) return(w) else return(r)
  }

#m = matrix(c(10,5,15,5,20,10,0,5,10,10,20,15,0,5,60,25),nrow=4,byrow=TRUE)
#dimnames(m) <- list(LETTERS[1:nrow(m)],as.roman(1:ncol(m)))
#m
#EU(m, p=c( .5, .25, .15, .1))
#maximin(m)
#hurwicz(m,1) # same as maximin
#
#m = matrix(c(10,5,15,5,20,10,0,5,10,10,20,15,0,5,60,25,21,8,2,11),nrow=5,byrow=TRUE)
#dimnames(m) <- list(LETTERS[1:nrow(m)],as.roman(1:ncol(m)))
#m
#EU(m, p=c( .5, .25, .15, .1))
#maximin(m)
#hurwicz(m) 
#chris_rock(m)
#minimaxregret(m)


##########################################################################
# Sensitivity analysis and tornado plots                                         [see marking.r for a more elaborated version]
##########################################################################

# The SA function offers just one approach to sensitivity analysis.  It is 1-var(pinched)/var(baseline) where 
# baseline is a distribution taking account of the variabilities of all inputs, and pinched is the distribution 
# for which one of the inputs has been set to its mean so that it contributes no variance to the calculation.
#
# https://en.wikipedia.org/wiki/Sensitivity_analysis#Variance-based_methods
# https://en.wikipedia.org/wiki/Sensitivity_analysis
# https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis

SA = function(x, f, tornado=TRUE, thick=0.60, col=1, border=NA, basecolor='gray98', skeleton=mean, wide=var) {
  k = length(x)
  t = paste0(deparse(f), collapse = " ")
  sf = deparse(f)[2]
  for (i in 1:k) sf = gsub(names(x[i]), paste('pa[[', i,']]',sep=''), sf, ignore.case = FALSE) 
  formals(f) = alist(pa=)
  body(f) <- str2expression(sf)
  ff = f(x)
  edf(ff)
  pl(ff); title(t)
  edf(ff,col=basecolor)
  m = mean(ff)
  base = var(ff)
  abline(v=m, lty='dotted', lwd=3)
  b = n = NULL
  for  (i in 1:k) {
    y = x
    y[[i]] = mean(y[[i]])
    n = c(n, names(x[i]))
    b = c(b, 1-var(f(y))/base) }
  bb = b[order(b)]
  rect(m-bb*base/5, ((1:k)-thick/2)/(k+1), m+bb*base/5, ((1:k)+thick/2)/(k+1),col=col,border=border)    
  text(left(ff),(1:k)/(k+1),n[order(b)],adj=c(0,0))
  names(b) = names(x)
  return(b)
  }
# Example
# R = normal(2, cv=0.05)          # hydraulic radius [m]
# S = normal(1, cv=0.10)          # slope of the energy line [%]
# n = normal(0.013, cv=0.10)    # coefficient of channel roughness
# f = function(R,S,n) R^(2/3) * S^(1/2) / n    
# s = SA(list(R=R, S=S, n=n), f)
# s = SA(list(R=R, S=S, n=n, dummy=0), f)


##########################################################################
# Distribution comparisons, validation, distance metrics and norms
##########################################################################

# See the gofTest function in the R package 'EnvStats' [https://search.r-project.org/CRAN/refmans/EnvStats/html/gofTest.html]
# See also the R package 'dgof' [https://cran.r-project.org/web/packages/dgof/dgof.pdf]

# see also the cor( ) function


#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#andersondarling2 = function(x,y,pl=FALSE,...) { # the A-D distance between two distributions
#  nx = length(x@x)
#  ny = length(y@x)
#  o = order(c(x@x,y@x))
#  r = c(rep(1,nx),rep(0,ny))
#  cx = cumsum(r[o])/nx
#  cy = cumsum((1-r)[o])/ny
#  okay = (cy != 0) & (cy != 1)
#  cx = cx[okay]
#  cy = cy[okay]
#  return(sum((cx - cy)^2 / (cy * (1-cy)))/1000)
#  } 
#
#andersondarling = function(A,B,pl=FALSE,...) { # the A-D distance between two distributions
#  ra = range(A@x)
#  rb = range(B@x)
#  Fn = F = NULL
#  for (i in seq(min(ra[[1]],rb[[1]]), max(ra[[2]],rb[[2]]), length.out=1000)) { Fn = c(Fn,prob(A,i)); F = c(F,prob(B,i)) }
#  okay = (F != 0) & (F != 1)
#  F = F[okay]
#  Fn = Fn[okay]
#  return(sum((Fn - F)^2 / (F * (1-F))))
#  } 
#  
#AD1 = AD2 = NULL
#for (i in 1:100) {
#  x = N(5,1)
#  y = U(3,9)
#  AD1 = c(AD1, andersondarling(x,y))
#  AD2 = c(AD2, andersondarling2(x,y))
#  }
#  
#pl(AD1,AD2)
#red(AD1)
#blue(AD2)
#
#
#  
#rbyc(5,5)  ############################
#mm = rep(c(100,250,500,1000,10000),5)
#MM = rep(5+c(-2,-1,0,1,2),rep(5,5))
#MC$many = 20000;   a  = N(5,1)
#ME = D = DD = NULL
#for (i in 1:25) {
#  m = mm[[i]]
#                                                                m = 20000
#									     M = MM[[i]]
#  MC$many = m;
#									     #b = N(6,1)
#									     #b  = N(M,1)
#									     me = runif(1,1,9)
#									     #b  = N(me,1)
#									     b = poisson(me)
#  pl(3,10)
#  edf(a)
#  blue(b,lwd=1)
#  d = andersondarling(a,b)
#  dd = andersondarling2(a,b)
#  D = c(D,d)
#  DD = c(DD,dd)
#  ME = c(ME,me)
#  title(paste('many =',m,'   A2 =',signif(d,3),signif(dd)))
#  }
#plot(abs(ME-5),log(D)) 
#plot(abs(ME-5),log(DD)) 
#
#cor(log(D),log(DD))
#cor(abs(ME-5),log(D)) 
#cor(abs(ME-5),log(DD)) 
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#ad = function(x,
#
#X = sort(histogram(x)@x)
#Y = sort(histogram(y)@x)
#
#
#
#
#
## ci doesn't need a comma between x and ... ??
#
#chisq.gof = function(x,y,df=25)  { # x is the observed; y is the expected
#  partition = c(min(y),cut.mc(y, seq(0,1,length.out=df+1)))  # seq(r[[1]], r[[2]], length.out=df+1)
#  s = 0
#  squeal = NULL
#  for (i in 1:df) {
#    o = sum(((partition[[i]] <= x@x) & (x@x < partition[[i+1]])) + 0)
#    e = sum(((partition[[i]] <= y@x) & (y@x < partition[[i+1]])) + 0)
#    if (e < 5) squeal = c(squeal,e)
#    s = s + (o - e)^2 / e
#    }
#  if (!is.null(squeal)) cat('chisq.gof Ei < 5\n',squeal,'\n')
#  return(s)
#  }
#
#O = c(81,50,27)
#res = chisq.test(O)   
#res$statistic
#
#E = rep(sum(O) / 3, 3)
#res = chisq.test(O, E)   
#res$statistic
#
#x = histogram(O)
#y = histogram(E)
#chisq.gof(x,y,df=1)  
#chisq.gof(x,y,df=2)  
#chisq.gof(x,y,df=3)  
#chisq.gof(x,y,df=4)  
#
#
#
#
#
#
#MC$many = 20000
#x = N(5,1)
#y = N(5,1)
#
#chisq.gof(x,y,df=15)  
#chisq.gof(x,y,df=25)  
#chisq.gof(x,y,df=50)  
#chisq.gof(x,y,df=200)  
#chisq.gof(x,y,df=1000)  
#
#
#chisq.gof(x,y,df=15) /15
#chisq.gof(x,y,df=25)  /25
#chisq.gof(x,y,df=50)  /50
#chisq.gof(x,y,df=200)  /200
#chisq.gof(x,y,df=1000)  /1000
#
#
#
#
#




# the areametric is NOT the same as the Wasserstein metric (they generalize to multiple dimensions differently)
areametric = function(x,y) sum(abs(sort(x@x) - sort(y@x))) / length(x@x)

smirnov = function(x,y,pl=FALSE,...) { # the Smirnov distance between two distributions, i.e., the maximal distance between their cdfs
  nx = length(x@x)
  ny = length(y@x)
  o = order(c(x@x,y@x))
  r = c(rep(1,nx),rep(0,ny))
  cx = cumsum(r[o])/nx
  cy = cumsum((1-r)[o])/ny
  if (pl) w = which.max(abs(cx-cy)); lines(rep(c(x@x,y@x)[o][[w]],2), c(cx[[w]], cy[[w]]), ...) # plots a line segment of length equal to the maximum distance, assuming x and y have already been plotted
  return(max(abs(cx-cy)))
  } 
#rbyc(5,5)  ############################
#mm = rep(c(100,250,500,1000,10000),5)
#MC$many = 20000;   a  = N(5,1)
#for (i in 1:25) {
#  m = mm[[i]]
#  MC$many = m;   b  = N(6,1)
#  pl(3,10)
#  edf(a)
#  blue(b,lwd=1)
#  d = smirnov(a,b,TRUE,col=2,lwd=2)
#  title(paste('many =',m,'   Dmax =',d))
#  }
#rbyc(5,5)  ############################
#mm = rep(c(100,250,500,1000,10000),5)
#MC$many = 20000;   a  = N(5,1)
#for (i in 1:25) {
#  m = mm[[i]]
#  MC$many = m;   b  = U(5,10)
#  pl(3,10)
#  edf(a)
#  lines(c(5,10),c(0,1),col='grey80')
#  blue(b,lwd=1)
#  d = smirnov(a,b,TRUE,col=2,lwd=2)
#  title(paste('many =',m,'   Dmax =',d))
#  }
#rbyc(5,5)  ############################
#mm = rep(c(100,250,500,1000,10000),5)
#MC$many = 20000;   a  = N(5,1)
#for (i in 1:25) {
#  m = mm[[i]]
#  MC$many = m;   b  = poisson(5)
#  pl(0,10)
#  edf(a)
#  blue(b,lwd=1)
#  d = smirnov(a,b,TRUE,col=2,lwd=2)
#  title(paste('many =',m,'Dmax =',d))
#  }
  
dom = function(b, c, tol=0.01)   { # does b stochastically dominate c?  reduce tol to make the comparison stricter, if neither dominates the other, they cross; if each does dominate the other, they're equal
  d = sort(mc(b)@x) - sort(mc(c)@x)
  if (tol==0) return(all(0<=d))
  return(-2*(sum(d[d<0])/(mean(b)+mean(c))) < tol * reps(b)) }
#
#################################################
# rbyc(5,4)
# a = N(5,1)
# Sdom = function(b,c,tol=0.01) {pl(b,c); blue(b); red(c); t = dom(b,c,tol); title(paste(t,' ... ', dom(c,b,tol))); t}
# # these will always be recognized as the same (TRUE ... TRUE)
# Sdom(a, a)
# Sdom(a, samedistribution(a))
# Sdom(a + 0.001, a)
# Sdom(N(5,1), a)
# # these are always recognized as crossing (FALSE ... FALSE)
# Sdom(uniform(1,9),a)
# Sdom(round(a),a)
# Sdom(N(5,0.85),a)
# Sdom(5,a)
# # these are almost always recognized as red dominating blue (FALSE ... TRUE)
# Sdom(N(4.6,1), a)
# Sdom(N(4.7,1), a)
# Sdom(N(4.8,1), a)
# Sdom(N(4.9,1), a)
# # these will usually be recognized as the same (TRUE ... TRUE) 
# Sdom(N(5.0,1), a); title('=')
# Sdom(N(5.0,1), a); title('=')
# Sdom(N(5.0,1), a); title('=')
# Sdom(N(5.0,1), a); title('=')
# # these are almost always recognized as blue dominating red (TRUE ... FALSE)
# Sdom(N(5.1,1), a)
# Sdom(N(5.2,1), a)
# Sdom(N(5.3,1), a)
# Sdom(N(5.4,1), a)
# ################################################
# # the same distribution should stochastically dominate itself
# rbyc(5,4)
# tol = 0.01      # almost all of these cases are  recognized as the same distribution (TRUE ... TRUE)
# m = runif(1,0,10)
# s = runif(1,0,10)
# a = N(m,s)
# for (kk in 1:20) Sdom(N(m,s), a, tol)
# c(m,s)
# ################################################
# rbyc(5,4)
# tol = 0.001; for (kk in 1:20) Sdom(N(m,s), a, tol)    # about half of the cases are recognized as the same
# ################################################
# rbyc(5,4)
# tol = 0.01    # almost all of the cases are recognized as the same
# m = runif(1,0,10)
# a = poisson(m)
# for (kk in 1:20) Sdom(poisson(m), a, tol)
# c(m)
# ################################################
# rbyc(5,4)
# tol = 0.005;  for (kk in 1:20) Sdom(poisson(m), a, tol)    # about half if the cases are recognized as the same
# ################################################


##########################################################################
# Replicate and sequence operations 
##########################################################################

# not to be confused with the reps( ) function

rep.mc <- function(x, ...) {
  r = NULL
  for (i in 1:list(...)[[1]]) r = c(r,x)
  return(r) 
  }


##########################################################################
# Exporting, writing, printing, plotting, summary, and str functions 
##########################################################################

summary.mc <- function (object, ...) {
  if (!is.mc(object)) stop('Object is not an MC distribution')
  ans <- list(mean='', variance='', sd='', iqr='', iqr.width=0, range='', left=0, 
   percentile.01='', percentile.05='', percentile.25='', median='',
   percentile.75='', percentile.95='', percentile.99='', right=0, replications=0)
  ans$mean              <- mean.mc(object)
  ans$variance          <- var.mc(object)
  ans$sd                   <- sd.mc(object)
  ans$iqr                  <- iqr.mc(object)
  #ans$iqr.width       <- width.interval(iqr.mc(object))
  ans$range              <- range.mc(object)
  ans$left                  <- left.mc(object)
  ans$percentile.01    <- cut.mc(object, 0.01)
  ans$percentile.05    <- cut.mc(object, 0.05)
  ans$percentile.25    <- cut.mc(object, 0.25)
  ans$median            <- cut.mc(object, 0.50)
  ans$percentile.75    <- cut.mc(object, 0.75)
  ans$percentile.95    <- cut.mc(object, 0.95)
  ans$percentile.99    <- cut.mc(object, 0.99)
  ans$right                <- right.mc(object)
  ans$replications      <- object@n
  class(ans)               <- "summary.mc"
  ans
  }

# specify other args by calling this DIRECTLY. as "print.summary.mc(summary.mc(x), digits=4)"

print.summary.mc <- function(x, ...) {  
  replaceempty <- function(s,r) if (s=='') r else s
  cat('\nMonte Carlo distribution summary',
      '\n  Mean: ',format(x$mean,...),
      '\n  Variance: ',format(x$variance,...),
      '\n  Std Deviation: ',format(x$sd,...),
      '\n  Width of interquartile range: ',format(diff(x$iqr),...),
      '\n  Width of overall range: ',format(diff(x$range),...),
      '\n  Order statistics',
      '\n     Left (min) value: ',format(x$left,...),
      '\n     1st percentile: ', format(x$percentile.01,...),
      '\n     5th percentile: ', format(x$percentile.05,...),
      '\n     25th percentile: ', format(x$percentile.25,...),
      '\n     Median (50th%ile): ', format(x$median,...),
      '\n     75th percentile: ', format(x$percentile.75,...),
      '\n     95th percentile: ', format(x$percentile.95,...),
      '\n     99th percentile: ', format(x$percentile.99,...),
      '\n     Right (max) value: ', format(x$right,...),
      '\n  Replications: ', format(x$replications,...),
      '\n',sep = '')
  invisible(x)
  }

as.character.mc <- function(x, ...) paste(sep='','MC (min=',min.mc(x),', median=',median.mc(x),', mean=',mean.mc(x),', max=',max.mc(x),')')

print.mc <- function(x, ...) cat(as.character.mc(x, ...),'\n')

show.mc <- function(x, ...) {
  try(plot.mc(x, ...))
  print.mc(x, ...)
  }
#show.mc <- function(x, xlab=deparse(substitute(x)), new=TRUE, ...) {
#  try(edf.mc(x, xlab=xlab, new=TRUE, ...))
#  print.mc(x, ...)
#  }
 
quiet <- setMethod("show", "mc", function(object)show.mc(object))
 
updown <- function(cumulative, x) if (cumulative) x else 1-x

cumulative <- function(yes=TRUE) MC$cumulative <- yes

#plot.mc <- function(s, xlab='', ylab='Cumulative probability', cumulative=MC$cumulative, ...) {
#  if (!cumulative) ylab = 'Exceedance probability'
#  plot(c(min(s@x),sort(s@x)),updown(cumulative,seq(0,1,length.out=1+length(s@x))),type='s', xlab=xlab,ylab=ylab,...)
#  bringToTop(-1) # return focus to R console
#  }

# if first two args are both mc, bivariate plot is made, else the edf is displayed using the second arg, if provided, as the label 
plot.mc <- function(x, t=NULL, xlab=t, ylab='Cumulative probability', cumulative=MC$cumulative, ...) {
  if (is.mc(t)) {
        xlab <- deparse(substitute(x))
        ylab <- deparse(substitute(t))
        plot(x@x,t@x,xlab=xlab,ylab=ylab, ...)
        } else { 
        if (missing(xlab)) xlab <- deparse(substitute(x))
        if (!cumulative) ylab = 'Exceedance probability'
        plot(c(min(x@x),sort(x@x)),updown(cumulative,seq(0,1,length.out=1+length(x@x))),type='s', xlab=xlab,ylab=ylab,...)
        }
  if (.Platform$OS.type=='windows') if(!RStudio) bringToTop(-1) # return focus to R console
  }
  
edf.numeric <- function(x, new=dev.cur()==1, xlab='', ylab='Cumulative probability', xlim=range(x), ylim=c(0,1), ...) {
  if (new) {plot(NULL,NULL,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim, ...); if (.Platform$OS.type=='windows') if(!RStudio) bringToTop(-1)}
  lines(c(min(x),sort(x)),seq(0,1,length.out=1+length(x)),type='s', ...)
  }
  
edf.mc <- function(b,xlab=deparse(substitute(b)),new=dev.cur()==1,...) edf.numeric(b@x,xlab=xlab,new=new,...)

range.list <- function(c, na.rm = FALSE) {r = range(c[[1]]); for (i in 2:length(c)) r = range(c(r,range(c[[i]]))); r}
#if(!isGeneric("range")) quiet <- setGeneric("range", function(x, ...) standardGeneric("range"))
#quiet <- setMethod('range', 'list', function(x, ...) range.list(x,...))

RANGE <- function(x) diff(range(x))

edf.list <- function(c,new=dev.cur()==1,xlim=range.list(c),...) for (i in 1:length(c)) edf(c[[i]], xlim=xlim, ...)

edf <- function(...) UseMethod('edf')
#if(!isGeneric("edf")) quiet <- setGeneric("edf", function(x, ...) standardGeneric("edf"))
#quiet <- setMethod('edf', 'list', function(x, ...) edf.list(x,...))
#quiet <- setMethod('edf', 'mc', function(x, ...) edf(x,...))
#quiet <- setMethod('edf', 'numeric', function(x, ...) edf.numeric(x,...))

lines.mc <- points.mc <- function(s,xlab='',cumulative=MC$cumulative, ...) {
  lines(c(min(s@x),sort(s@x)),updown(cumulative,seq(0,1,length.out=1+length(s@x))),type='s',...)
  if (.Platform$OS.type=='windows') if(!RStudio) bringToTop(-1) # return focus to R console
  }

plot.mc.scale <- function(min, max, name='',cumulative=MC$cumulative, col=NULL, ...) {
  if (cumulative) ylab = 'Cumulative probability' else ylab = 'Exceedance probability'
  plot(c(min,max),c(0,1),xlab=name, ylab=ylab, col='white', ...)
  }

plot.mclist <- function(A,...) { # see also pledf( ) which also handles number arrays
  rm <-  Inf; for (i in 1:length(A)) {r <-  left.mc(A[[i]]); if (is.finite(r)) rm <- base::min(r,rm) }
  rM <- -Inf; for (i in 1:length(A)) {r <- right.mc(A[[i]]); if (is.finite(r)) rM <- base::max(r,rM) }
  if ("xlim" %in% attributes(list(...))$names) plot.mc.scale(list(...)$xlim[1],list(...)$xlim[2],...) else plot.mc.scale(rm, rM, ...)    
  for (i in 1:length(A)) lines(A[[i]], ...)
  }

#pl <- function(...) plot(NULL,ylim=c(0,1),xlim=range(...),xlab='',ylab='Cumulative probability')
#pl <- function(...,xlab='',ylab='Cumulative probability') plot(NULL,ylim=c(0,1),xlim=range(...),xlab=xlab,ylab=ylab)
#pl = function(..., xlab='',ylab='Cumulative probability',ylim=c(0,1)) plot(NULL,ylim=ylim,xlim=range(...),xlab=xlab,ylab=ylab)
pl = function(..., xlab='',ylab='Cumulative probability',ylim=c(0,1)) plot(NULL,ylim=ylim,xlim=range.mc(...),xlab=xlab,ylab=ylab) # see range.mc
pledf = function(...) {pl(...); elts=list(...); if (mode(elts[[1]])=='list') elts=elts[[1]]; for (e in elts) edf(e)} #display several objects at once (to display one MC object, just enter its name)
rbyc = function(r=1,c=1) par(mfrow=c(r,c))
cbyr = function(c=1,r=1) par(mfcol=c(r,c))

black      = bla = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='black',...)      edf(x,xlab=xlab,lwd=lwd,col=col,...)
red        = red = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='red',...)        edf(x,xlab=xlab,lwd=lwd,col=col,...)
blue       = blu = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='blue',...)       edf(x,xlab=xlab,lwd=lwd,col=col,...)
green      = gre = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='green3',...)     edf(x,xlab=xlab,lwd=lwd,col=col,...)
khaki      = kha = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='khaki2',...)     edf(x,xlab=xlab,lwd=lwd,col=col,...)
navy       = nav = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='navy',...)       edf(x,xlab=xlab,lwd=lwd,col=col,...)
purple     = pur = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='purple',...)     edf(x,xlab=xlab,lwd=lwd,col=col,...)
brown      = bro = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='brown4',...)     edf(x,xlab=xlab,lwd=lwd,col=col,...)
sienna     = sie = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='sienna',...)     edf(x,xlab=xlab,lwd=lwd,col=col,...)
rust       = rus = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='sienna3',...)     edf(x,xlab=xlab,lwd=lwd,col=col,...)
olivedrab  = oli = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='olivedrab4',...) edf(x,xlab=xlab,lwd=lwd,col=col,...)
gray       = gra = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='gray',...)       edf(x,xlab=xlab,lwd=lwd,col=col,...)
orange     = ora = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='darkorange',...) edf(x,xlab=xlab,lwd=lwd,col=col,...)
pink       = pin = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='deeppink',...)   edf(x,xlab=xlab,lwd=lwd,col=col,...)
salmon     = sal = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='pink1',...)   edf(x,xlab=xlab,lwd=lwd,col=col,...)
magenta    = mag = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='magenta',...)   edf(x,xlab=xlab,lwd=lwd,col=col,...)
cyan       = cya = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='cyan',...)       edf(x,xlab=xlab,lwd=lwd,col=col,...)
chartreuse = cha = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='chartreuse',...) edf(x,xlab=xlab,lwd=lwd,col=col,...)
yellow     = yel = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='gold',...)       edf(x,xlab=xlab,lwd=lwd,col=col,...)
tan        = ecr = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='tan',...)        edf(x,xlab=xlab,lwd=lwd,col=col,...)
brick      = bri = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='firebrick3',...) edf(x,xlab=xlab,lwd=lwd,col=col,...)
white      = whi = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col='white',...) edf(x,xlab=xlab,lwd=lwd,col=col,...)
teal       = teal = function(x,xlab=deparse(substitute(x)),lwd=MC$lwd,col=rgb(0,128/255,128/255),...)       edf(x,xlab=xlab,lwd=lwd,col=col,...)

units = function(s) return(1) # the R library doesn't support units or units checking like RAMAS Risk Calc does

if ("MC$scott.used.to.be.quiet" %in% ls()) {quiet <- MC$scott.used.to.be.quiet; rm("MC$scott.used.to.be.quiet")} else rm("quiet")
cat(':sra> library loaded\n')

##########################################################################
# End of Monte Carlo Simulation S4 Library for the R Language                
##########################################################################

##########################################################################
# Some further programmer notes
##########################################################################

## The expression pmin(5, N(5,1) was not working, and the solution for the analogous 
## problem that we used for the infix operators didn't seem to be available.
## This was the error:
#
#pmin(N(5,1),N(6,2))                        # works with R's ordinary dispatch rules
#pmin(N(5,1),5)                             # works with R's ordinary dispatch rules
#pmin(5, N(5,1))                            # FAILS under R's ordinary dispatch rules
#pmin(c(2,3,4,16,17,18),c(12,13,14,6,7,8))  # must continue to work
#
## This is what turned out to be the fix:
#
#pmin.mc = function (..., na.rm = FALSE) {  
#cat('yo')                 ################### check that control gets here
#  elts <- makemc(...)
#  m <- elts[[1]]
#  for (each in elts[-1]) m <- conv.mc(m, each, 'pmin')
#  m
#  }
#
#pmin <- function(...,na.rm=FALSE) {UseMethod("pmin")}; 
#pmin.default <- function(..., na.rm = FALSE) base::pmin(...,na.rm=na.rm)
#pmin.numeric = function(...,na.rm=FALSE){e=list(...); if (class(e[[2]])=='mc') return(pmin.mc(mc(e[[1]]),e[[2]],na.rm=na.rm)) else return(base::pmin(...,na.rm=na.rm))}
#
## The problem disappeared:
#
#pmin(N(5,1),N(6,2))                        # works with R's ordinary dispatch rules
#pmin(N(5,1),5)                             # works with R's ordinary dispatch rules
#pmin(5, N(5,1))                            # NOW WORKS after .default and .numeric and .mc functions are defined
#pmin(c(2,3,4,16,17,18),c(12,13,14,6,7,8))  # must continue to work
#
## The following three approaches did not work
#
##infinite recursion
##pmin.numeric = function(...,na.rm=FALSE){e=list(...); if (class(e[[2]])=='mc') return(pmin.mc(mc(e[[1]]),e[[2]],na.rm=na.rm)) else UseMethod("pmin")}
#
##never gets to pmin.mc
##pmin <- function(...,na.rm=FALSE) {UseMethod("pmin")}; pmin.default <- function(..., na.rm = FALSE) base::pmin(...,na.rm=na.rm)
##pmin(5, N(5,1))
#
##infinite recursion
##pmin <- function(...,na.rm=FALSE) {UseMethod("pmin")}; pmin.default <- function(..., na.rm = FALSE) pmin.mc(...,na.rm=na.rm)
##pmin(5, N(5,1))

##########################################################################
# Handling implicit dependence
#
#a = N(5,2)
#pl(-5,25)
#
# shuffling
#cyan(a + samedistribution(a))
#red(selfsum(a,2))
#black(a + shuffle(a))
#
# preserves implicit dependence
#green(perfectconv.mc(a,a))
#blue(a * 2)
#gray(a+a)
 
