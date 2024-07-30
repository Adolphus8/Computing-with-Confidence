source("pba BETTER.r")

#######################################################################
readkeygraph <- function(prompt) {
getGraphicsEvent(prompt=prompt,onMouseDown=NULL,onMouseMove=NULL,onMouseUp=NULL,onKeybd=onKeybd,consolePrompt="[click on graph then follow top prompt to continue]")
Sys.sleep(0.01)
return(keyPressed)}

onKeybd <- function(key) keyPressed <<- key
wait = function(msg='') invisible(readkeygraph(paste(msg,'[press any key to continue]')))
rbyc = function(r=1,c=1) par(mfrow=c(r,c))
pl = function(..., xlab='',ylab='Cumulative probability',ylim=c(0,1)) plot(NULL,ylim=ylim,xlim=range(...),xlab=xlab,ylab=ylab)
sh = function(w, wi=NULL, m, pl=range(w)) {plot(w,xlim=pl); if (!missing(wi)) green(wi); blue(w); title(m)}
#######################################################################
samp = 200;

# Define the component C-boxes:
C1 = KN(23,24);
C2 = KN(23,24);
C3 = KN(14,17);
C4 = KN(14,17);
C5 = KN(12,12);

# Plot the C-boxes of the components:
rbyc(2,3)
sh(C1,,'Component 1',c(0,1))
sh(C2,,'Component 2',c(0,1))
sh(C3,,'Component 3',c(0,1))
sh(C4,,'Component 4',c(0,1))
sh(C5,,'Component 5',c(0,1))

# Define the Event C-boxes under Independence via Boolean logic (Path-set):
phi1 = andI(C1,C3);
phi2 = andI(C2,C4);
phi3 = andI(C1,andI(C4,C5));
phi4 = andI(C2,andI(C3,C5));
R_S1_ps_Boolean = orI(orI(phi1,phi2),orI(phi3,phi4));

source("sra.r")

#######################################################################
many = 10000
constant <- function(b) if (length(b)==1) TRUE else FALSE
precise <- function(b) if (length(b)==many) TRUE else FALSE
leftside <- function(b) if (precise(b)) return(b) else return(b[1:many])
rightside <- function(b) if (precise(b)) return(b) else return(b[(many+1):(2*many)])
sampleindices = function() round(runif(many)*(many-1) + 1)
env <- function(x,y) if ((precise(x) && precise(y))) c(x,y) else stop('env error') # only works for precise distributional inputs
beta <- function(v,w) if ((v==0) && (w==0)) env(rep(0,many),rep(1,many)) else if (v==0) rep(0,many) else if (w==0) rep(1,many) else sort(rbeta(many, v, w))
pairsides <- function(b) {i = sampleindices(); return(env(leftside(b)[i],rightside(b)[i]))}

kn <- function(k,n) return(pairsides(env(beta(k, n-k+1), beta(k+1, n-k)))) 
orI <- function(x,y) return(1-(1-x)*(1-y))
andI <- function(x,y) return(x*y)
#######################################################################

# Define the State function of the Bridge structure:
bb = function(a,b,c,d,e) a*c+b*d+a*d*e+b*c*e-a*b*c*e-a*c*d*e-a*b*d*e-b*c*d*e-a*b*c*d+2*a*b*c*d*e

bb_boolean = function(a,b,c,d,e) {
phi1 = andI(a,c);
phi2 = andI(b,d);
phi3 = andI(a,andI(d,e));
phi4 = andI(b,andI(c,e));
return(orI(orI(phi1,phi2),orI(phi3,phi4)));
}

bridgeI = function(a,b,c,d,e,bb) {
  L = leftside
  R = rightside
  return(env(bb(L(a),L(b),L(c),L(d),L(e)), bb(R(a),R(b),R(c),R(d),R(e))))
}

plotbox <- function(b,new=TRUE,col='blue',lwd=2,xlim=range(b[is.finite(b)]),ylim=c(0,1),xlab='',ylab='Cumulative probability',...) {

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

C1 = kn(23,24);
C2 = kn(23,24);
C3 = kn(14,17);
C4 = kn(14,17);
C5 = kn(12,12);

R_S1_ps_State = bridgeI(C1,C2,C3,C4,C5,bb);
R_S1_ps_Boolean = bridgeI(C1,C2,C3,C4,C5,bb_boolean);


# Plot the C-boxes of the system reliability for S1:
plotbox(R_S1_ps_Boolean); title('S1 Reliability - Boolean Function')
plotbox(R_S1_ps_State,FALSE,col='green'); title('S1 Reliability - State Function')
