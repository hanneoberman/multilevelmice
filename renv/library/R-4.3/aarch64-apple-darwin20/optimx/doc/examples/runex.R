# runex.R  -- Run selected scripts in the inst/doc/examples/
#    directory of package optimx
# J C Nash 2023-6-21
currdir <- getwd()
cat("Current directory is ",currdir,"\n")
cat("If this is the inst/doc/examples/ directory of package 'optimx'\n")
cat("then Enter, else provide the path to that directory\n");
newdir<-readline("Path to dir:")

cat(length(newdir),"\n")

if (length(newdir) > 1){
   setwd(newdir)
}

myrun<-function(cmdstr){
   cat("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")
   cat("About to run script file ",cmdstr,"\n")
   system(cmdstr)
   cat(">>>>>> end of ",cmdstr," >>>>>>>>>>>>>\n")
   tmp <- readline("Continue?")
   if (tmp == "STOP") stop("DONE!")
   return(0)
}


myrun("Rscript 3Rosen.R")
myrun("Rscript argclash.R")
myrun("Rscript axsearch-Ex.R")
myrun("Rscript brown_test.R")
myrun("Rscript broydt_test.R")
myrun("Rscript chenlog_test.R")
myrun("Rscript chen_test.R")
myrun("Rscript cyq_test.R")
myrun("Rscript dropnmbk.R")
myrun("Rscript froth_test.R")
myrun("Rscript genrose_test.R")
myrun("Rscript hessapptest.R")
myrun("Rscript hessian-used.R")
myrun("Rscript hobbs.R")
myrun("Rscript jonesrun.R")
myrun("Rscript maxtestJN.R")
myrun("Rscript ncgtests.R")
myrun("Rscript onepar_test.R")
myrun("Rscript optimrgrapprox.R")
myrun("Rscript ox.R")
myrun("Rscript poissmix_test.R")
myrun("Rscript rosbkext_test.R")
myrun("Rscript rosenbrock.R")
myrun("Rscript sc2_test.R")
myrun("Rscript scalechk-Ex.R")
myrun("Rscript simplefuntst.R")
myrun("Rscript snewtonbtest.R")
myrun("Rscript specctrlhobbs.R")
myrun("Rscript ssqtest.R")
myrun("Rscript trig1507.R")
myrun("Rscript trystarttests.R")
myrun("Rscript valley_test.R")
myrun("Rscript vmmix_test.R")
myrun("Rscript woodtest.R")
