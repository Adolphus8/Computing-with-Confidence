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

# Define the component C-boxes:
T = KN(0,2000);
K2 = KN(3,500);
S = KN(1,150);
S1 = KN(0,460);
K1 = KN(7, 10000);
R = KN(0,380);

# Define the Event C-boxes under uncertain dependencies:
e5f = or(K1,R);
e4f = or(e5f,S1);
e3f = and(e4f, S);
e2f = or(e3f, K2);
e1f = or(e2f, T);

# Define the Event C-boxes under independence:
e5i = orI(K1,R);
e4i = orI(e5i,S1);
e3i = andI(e4i, S);
e2i = orI(e3i, K2);
e1i = orI(e2i, T);

# Plot the C-boxes of the components:
rbyc(2,3)
sh(T,,'Pressure tank T',c(0,2*right(mean(T))))
sh(K2,,'Relay K2',c(0,2*right(mean(K2))))
sh(S,,'Pressure switch S',c(0,2*right(mean(S))))
sh(S1,,'On-off switch S1',c(0,2*right(mean(S1))))
sh(K1,,'Relay K1',c(0,2*right(mean(K1))))
sh(R,,'Time relay R',c(0,2*right(mean(R))))

wait()

# Plot the C-boxes of the Events
sh(e5f,e5i,'E5',c(0,2*right(mean(e5f))))
sh(e4f,e4i,'E4',c(0,2*right(mean(e4f))))
sh(e3f,e3i,'E3',c(0,2*right(mean(e3f))))
sh(e2f,e2i,'E2',c(0,2*right(mean(e2f))))
sh(e1f,e1i,'E1',c(0,2*right(mean(e1f))))

wait()

## Results of the 95% confidence intervals:
conf_lb = 0.025; conf_ub = 1 - conf_lb;

# Under independence:
interval_e1i = interval(cut(e1i, conf_lb), cut(e1i, conf_ub))
interval_e2i = interval(cut(e2i, conf_lb), cut(e2i, conf_ub))
interval_e3i = interval(cut(e3i, conf_lb), cut(e3i, conf_ub))
interval_e4i = interval(cut(e4i, conf_lb), cut(e4i, conf_ub))
interval_e5i = interval(cut(e5i, conf_lb), cut(e5i, conf_ub))

# Under frechet:
interval_e1f = interval(cut(e1f, conf_lb), cut(e1f, conf_ub))
interval_e2f = interval(cut(e2f, conf_lb), cut(e2f, conf_ub))
interval_e3f = interval(cut(e3f, conf_lb), cut(e3f, conf_ub))
interval_e4f = interval(cut(e4f, conf_lb), cut(e4f, conf_ub))
interval_e5f = interval(cut(e5f, conf_lb), cut(e5f, conf_ub))
 