
# Adolphus:
#
# You probably thought I was trying to drive you insane.  
#
# I discovered (in those old comments) the reason why the Boolean expresssion and its arithmetic expansion were different and why the left-right-Monte-Carlo approach does not work for the Boolean expression.  That is satisfying I think.
#
# However, I thought I'd also found the solution to the numerical discrepancy, but I was only at a local peak where, because I made errors in both codes in the same way, I was getting agreement.  A case of two wrongs seeming to make a right.
#
# Now, however, I think I've got good agreement, and for the right reasons.  This involved restoring the use of the sampleindices() function in the OG code.  It has the effect of randomly shuffling the deviates to correctly model independence.  
#
# Use File/Source to first load sra.r and then this file


bothmany = 100000
many = bothmany

KNL = function(k,n) rbeta(many, k, n - k + 1)
KNR = function(k,n) rbeta(many, k + 1, n - k)

C1L = KNL(23,24)
C2L = KNL(23,24) # do NOT use R's multiple assignment (C1L = C2L = KNL...) as that induces perfect dependence between them
C3L = KNL(14,17)
C4L = KNL(14,17)
C5L = KNL(12,12)

C1R = KNR(23,24)
C2R = KNR(23,24)
C3R = KNR(14,17)
C4R = KNR(14,17)
C5R = KNR(12,12)

par(mfrow=c(2,3))
pl(0,1); red(C1L); blue(C1R)
pl(0,1); red(C2L); blue(C2R)
pl(0,1); red(C3L); blue(C3R)
pl(0,1); red(C4L); blue(C4R)
pl(0,1); red(C5L); blue(C5R)

L = 1 - ((1 - C1L * C3L) * (1 - C2L*C4L) * (1 - C1L*C4L*C5L) * (1-C2L*C3L*C5L))
R = 1 - ((1 - C1R * C3R) * (1 - C2R*C4R) * (1 - C1R*C4R*C5R) * (1-C2R*C3R*C5R))

pl(0,1)
red(L)
blue(R)

x1 = C1L
x2 = C2L
x3 = C3L
x4 = C4L
x5 = C5L

xx = function(x1,x2,x3,x4,x5) (x1 * x3) + 
(x2 * x4) + 
(x1 * x4 * x5) +
(x2 * x3 * x5) - 
(x1 * x2 * x3 * x5) - 
(x1 * x3 *x4 * x5) - 
(x1* x2 * x4 * x5) - 
(x2 * x3 * x4 * x5) - 
(x1 * x2 * x3 * x4) +
(2 * x1 * x2 * x3 * x4 * x5)

L = xx(x1,x2,x3,x4,x5)

x1 = C1R
x2 = C2R
x3 = C3R
x4 = C4R
x5 = C5R

R = xx(x1,x2,x3,x4,x5)

pl(0,1)
red(L,lwd=1)
blue(R,lwd=1)
range(ci(L), ci(R))

# Reliability uncertainty with confidence boxes	
# https://sites.google.com/site/reliabilityuncertainty/code/8-bridge-relindep
# 8 Bridge Rel(indep)

###############################################################################
# Compute the reliability of a bridge with structure whose Boolean expression necessarily has repeated variables
# rm(list=ls())     # wipe everything from memory
# The bridge system is
#
#                +---a----+----c---+
#                |        |        | 
#          ------|        e        |-------
#                |        |        |
#                +---b----+----d---+ 
#
# where the reliabilities of the five (unrepairable) components are 
# given as independent uncertain numbers.  This calculation uses  
# a combination of monotonicity and Monte Carlo methods to  
# take account of independence assumptions and sidestep  
# the problem of the repeated variables.  The Boolean function is
# (a&c)|(b&d)|(a&e&d)|(b&e&c) which has irreducible repetitions
# of every variable.  Simplification to an arithmetic function yields
# ac+bd+ade+bce-abce-acde-abde-bcde-abcd+2abcde.
# a*c+b*d+a*d*e+b*c*e-a*b*c*e-a*c*d*e-a*b*d*e-b*c*d*e-a*b*c*d+2*a*b*c*d*e

# infrastructure
many = bothmany
constant <- function(b) if (length(b)==1) TRUE else FALSE
precise <- function(b) if (length(b)==many) TRUE else FALSE
leftside <- function(b) if (precise(b)) return(b) else return(b[1:many])
rightside <- function(b) if (precise(b)) return(b) else return(b[(many+1):(2*many)])
sampleindices = function() round(runif(many)*(many-1) + 1) # beta produces sorted deviates, so sampleindices restores independence in the pairsides function    
env <- function(x,y) if ((precise(x) && precise(y))) c(x,y) else stop('env error') # only works for precise distributional inputs
beta <- function(v,w) if ((v==0) && (w==0)) env(rep(0,many),rep(1,many)) else if (v==0) rep(0,many) else if (w==0) rep(1,many) else sort(rbeta(many, v, w))
pairsides <- function(b) {i = sampleindices(); return(env(leftside(b)[i],rightside(b)[i]))}
plotbox <- function(b,new=TRUE,col='blue',lwd=2,xlim=range(b[is.finite(b)]),ylim=c(0,1),xlab='',ylab='Probability',...) {
  edf <- function (x, col, lwd, ...) {
      n <- length(x)
      s <- sort(x)
      if (2000<n) {s = c(s[1],s[1:(n/5)*5],s[n]); n <- length(s);} # subsetting for the plot to every 5th value
      lines(c(s[[1]],s[[1]]),c(0,1/n),lwd=lwd,col=col,...)
      for (i in 2:n) lines(c(s[[i-1]],s[[i]],s[[i]]),c(i-1,i-1,i)/n,col=col,lwd=lwd,...)
      }
  b = ifelse(b==-Inf, xlim[1] - 10, b)
  b = ifelse(b==Inf, xlim[2] + 10, b)
  if (new) plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
  if (length(b) < many) edf(b,col,lwd) else
  edf(c(min(b),max(b),b[1:min(length(b),many)]),col,lwd)
  if (many < length(b)) edf(c(min(b),max(b),b[(many+1):length(b)]),col,lwd)
  }

  
# c-box for Bernoulli parameter:  x[i] ~ Bernoulli(parameter), x[i] is either 0 or 1, n = length(x), k = sum(x)
kn <- function(k,n) return(pairsides(env(beta(k, n-k+1), beta(k+1, n-k)))) 

bb = function(a,b,c,d,e) a*c+b*d+a*d*e+b*c*e-a*b*c*e-a*c*d*e-a*b*d*e-b*c*d*e-a*b*c*d+2*a*b*c*d*e

bridgeI = function(a,b,c,d,e) {
  L = leftside
  R = rightside
  return(env(bb(L(a),L(b),L(c),L(d),L(e)), bb(R(a),R(b),R(c),R(d),R(e))))
  }

a = kn(23,24)   
b = kn(23,24)  
c = kn(14,17) 
d = kn(14,17)
e = kn(12,12)
z = bridgeI(a,b,c,d,e)
# mean of the z structure?
c(mean(z[1:many]), mean(z[(many+1):(2*many)]))
par(mfrow=c(2,3))
plotbox(a); title('a'); green(C1L); green(C1R)
plotbox(b); title('b'); green(C2L); green(C2R)
plotbox(c); title('c'); green(C3L); green(C3R)
plotbox(d); title('d'); green(C4L); green(C4R)
plotbox(e); title('e'); green(C5L); green(C5R)
plotbox(z,lwd=1); title('bridge (assuming independence)'); green(L,lwd=1); green(R,lwd=1)

rbyc(); plotbox(z,lwd=1); title('bridge (assuming independence)'); green(L,lwd=1); green(R,lwd=1)
  
  
# three or four sig figs with many = 100000
   
mean(L)   #  0.9423821
mean(z[1:many])     #  0.9424204
 
mean(R)    #  0.9707185
mean(z[(many+1):(2*many)])  #  0.9702578
  
sd(L)   # 0.03253915
sd(z[1:many])     # 0.03285308

sd(R)    # 0.02158358
sd(z[(many+1):(2*many)])  # 0.0220665

ci(L)    #  0.8615944 0.9855997
ci(R)    #  0.9143150 0.9956702
range(ci(leftside(z)), ci(rightside(z)))   #  0.8596932 0.9956000
 





#A=runif(1); B=runif(1); C=runif(1); D=runif(1); E=runif(1)
#xx(A,B,C,D,E)
#bb(A,B,C,D,E)


cat('End of the calculation\nWe induce an error to stop reading the file\n')

++ # this error is intentional (so that you don't run the rest of the file without seeing what's displayed by the first part of this file)
  
  
  
  
  
  
# reproduce the interval calculation from <<citation>>
a = c(rep(0.7,many), rep(0.8,many))
b = c(rep(0.75,many), rep(0.9,many)) 
c = c(rep(0.9,many), rep(0.95,many)) 
d = 0.95
e = 0.82
zz = bridgeI(a,b,c,d,e)  #  Interval: [0.91556, 0.975327]

plotbox(a); title('a'); green(C1L); green(C1R)
plotbox(b); title('b'); green(C2L); green(C2R)
plotbox(c); title('c'); green(C3L); green(C3R)
plotbox(d); title('d'); green(C4L); green(C4R)
plotbox(e); title('e'); green(C5L); green(C5R)
plotbox(zz)
bb(.7,.75,.9,.95,.82)   # 91556
bb(.8,.9,.95,.95,.82)  # 0.975327

# reproduce output picture here <<is the result z correct?>>

###############################################################################

# The Boolean function (a&c)|(b&d)|(a&e&d)|(b&e&c) has repetitions, but it 
# cannot be evaluated using the same paired-sides approach that we use for
# the arithmetic function ac+bd+ade+bce-abce-acde-abde-bcde-abcd+2abcde.
# The plots below illustrate why.  The addends of the arithmetic expression 
# correspond to non-overlapping regions of the Venn diagram, so they are
# mutually exclusive. The disjunctands of the Boolean expression, however,
# are dependent of one another, but in a way that the Monte Carlo 
# calculation does not fix.

bb = function(a,b,c,d,e) a*c+b*d+a*d*e+b*c*e-a*b*c*e-a*c*d*e-a*b*d*e-b*c*d*e-a*b*c*d+2*a*b*c*d*e

orI <- function(x,y) return(1-(1-x)*(1-y))
andI <- function(x,y) return(x*y)

BB = function(a,b,c,d,e) orI(orI(andI(a,c), andI(b,d)), orI(andI(a, andI(e,d)), andI(b, andI(e,c))))

BB = function(a,b,c,d,e) orI(orI(a*c, b*d), orI(a*e*d, b*e*c))

BB2 = function(a,b,c,d,e) {
  Many = 1000
  n = length(a)
  w = rep(0, n)
  for (i in 1:Many) {
    A = runif(n) < a
    B = runif(n) < b
    C = runif(n) < c
    D = runif(n) < d
    E = runif(n) < e
    w = w + orI(orI(A*C, B*D), orI(A*E*D, B*E*C))
    }
  w / Many
}

bridgeI = function(a,b,c,d,e,bb) {
  L = leftside
  R = rightside
  return(env(bb(L(a),L(b),L(c),L(d),L(e)), bb(R(a),R(b),R(c),R(d),R(e))))
}

z = bridgeI(a,b,c,d,e,bb)
y = bridgeI(a,b,c,d,e,BB)
x = bridgeI(a,b,c,d,e,BB2)
par(mfrow=c(1,1))
plotbox(z)
plotbox(y, new = FALSE, col='gray')
plotbox(x, new = FALSE, col='red')
title('bridge (assuming independence)')