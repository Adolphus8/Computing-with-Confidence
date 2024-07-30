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
samp = 200;
p = KN(1,42);
v3 = KN(8,42);
ov3 = pbox(3.49e-5);
v4 = KN(8,42);
ov4 = pbox(3.49e-5);

# Define the Event C-boxes under uncertain dependencies:
e4f = or(v4,ov4);  cbox_vec = array(c(e4f@d, e4f@u),dim = c(samp,2)); write.csv(cbox_vec, "E4_cbox_frechet_TRIGA.csv")  
e3f = or(v3,ov3);  cbox_vec = array(c(e3f@d, e3f@u),dim = c(samp,2)); write.csv(cbox_vec, "E3_cbox_frechet_TRIGA.csv")  
e2f = or(e3f,e4f); cbox_vec = array(c(e2f@d, e2f@u),dim = c(samp,2)); write.csv(cbox_vec, "E2_cbox_frechet_TRIGA.csv")  
e1f = and(p, e2f); cbox_vec = array(c(e1f@d, e1f@u),dim = c(samp,2)); write.csv(cbox_vec, "E1_cbox_frechet_TRIGA.csv")  

# Define the Event C-boxes under independence:
e4i = orI(v4,ov4);  cbox_vec = array(c(e4i@d, e4i@u),dim = c(samp,2)); write.csv(cbox_vec, "E4_cbox_indep_TRIGA.csv") 
e3i = orI(v3,ov3);  cbox_vec = array(c(e3i@d, e3i@u),dim = c(samp,2)); write.csv(cbox_vec, "E3_cbox_indep_TRIGA.csv") 
e2i = orI(e3i,e4i); cbox_vec = array(c(e2i@d, e2i@u),dim = c(samp,2)); write.csv(cbox_vec, "E2_cbox_indep_TRIGA.csv") 
e1i = andI(p, e2i); cbox_vec = array(c(e1i@d, e1i@u),dim = c(samp,2)); write.csv(cbox_vec, "E1_cbox_indep_TRIGA.csv") 

# Plot the C-boxes of the components:
rbyc(2,3)
sh(p,,'Inlet pipeline',c(0, 0.35))
sh(v3,,'Valve V-3',c(0, 0.6))
sh(ov3,,'Operator of V-3',c(2e-5, 5e-5))
sh(v4,,'Valve V-4',c(0, 0.6))
sh(ov4,,'Operator of V-4',c(2e-5, 5e-5))

wait()

# Plot the C-boxes of the Events
rbyc(2,2)
sh(e1f,e1i,'E1',c(0,0.2))
sh(e2f,e2i,'E2',c(0,1))
sh(e3f,e3i,'E3',c(0,0.4))
sh(e4f,e4i,'E4',c(0,0.4))

wait()

## Results of the 95% confidence intervals:
conf_lb = 0.025; conf_ub = 1 - conf_lb;

# Under independence:
interval_e1i = interval(cut(e1i, conf_lb), cut(e1i, conf_ub))
interval_e2i = interval(cut(e2i, conf_lb), cut(e2i, conf_ub))
interval_e3i = interval(cut(e3i, conf_lb), cut(e3i, conf_ub))
interval_e4i = interval(cut(e4i, conf_lb), cut(e4i, conf_ub))

# Under frechet:
interval_e1f = interval(cut(e1f, conf_lb), cut(e1f, conf_ub))
interval_e2f = interval(cut(e2f, conf_lb), cut(e2f, conf_ub))
interval_e3f = interval(cut(e3f, conf_lb), cut(e3f, conf_ub))
interval_e4f = interval(cut(e4f, conf_lb), cut(e4f, conf_ub))
