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

# C-boxes of components:
C1 = KN(23,24) # Component 1
C2 = KN(14,17) # Component 2

rbyc(1,2)
C1; title('One')
C2; title('Two')

## Evaluate Series configuration: 
Ser = and(C1,C2)
Seri = andI(C1,C2)
Ser

green(Seri)
blue(Ser)

title('Series')

wait()

# Evaluate results:
conf_lb = 0.025; conf_ub = 1 - conf_lb;

interval_ser_i = interval(cut(Seri, conf_lb), cut(Seri, conf_ub)) # Under independence
interval_ser_f = interval(cut(Ser, conf_lb), cut(Ser, conf_ub))   # Under uncertain dependency

## Evaluate Parallel configuration:
Par = or(C1,C2)
Pari = orI(C1,C2)
Par

green(Pari)
blue(Par)

title('Parallel')

wait()

# Evaluate results:
conf_lb = 0.025; conf_ub = 1 - conf_lb;

interval_par_i = interval(cut(Pari, conf_lb), cut(Pari, conf_ub)) # Under independence
interval_par_f = interval(cut(Par, conf_lb), cut(Par, conf_ub))   # Under uncertain dependency
