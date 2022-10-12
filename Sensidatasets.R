##------SENSITIVITY ANALYSIS-----------------
#----------STEP 1: GETTING 10DATASETS FROM OUTDEOPTIM FOR 2015-2021 AND 2010-2021 FOR LOC_FERISSBURG-----------
#LOC_Ferrisburg_2010, is the  name of folder for LOC_ferrisburg 2010-2021 SWAT initialization
#LOC_Ferrisburg_2015, is the  name of folder for LOC_ferrisburg 2015-2021 SWAT initialization
dir.create("~/multibasinpaper")
setwd("~/multibasinpaper")
library(SWATmodel)
flowgage_id="04282650" #Little Otter Creek at Ferrisburg, VT.
flowgage=get_usgs_gage(flowgage_id,begin_date = "2010-01-01",end_date= "2022-01-01")
flowgage$flowdata$Qm3ps= flowgage$flowdata$flow/24/3600 #m3/s
######remove extreme datapoint
max(flowgage$flowdata$Qm3ps)
flowgage$flowdata=subset(flowgage$flowdata, flowgage$flowdata$Qm3ps<30)

##moving initialization folder to dev/shm directory for faster run
dir.create("/dev/shm/rojakaveh")
setwd("/dev/shm/rojakaveh")
dir.create("CrackFlowFer2010")
dir.create("CrackFlowFer2015")
file.copy(list.files("~/multibasinpaper/LOC_Ferrisburg_2010/",full.names = TRUE,recursive = TRUE),"/dev/shm/rojakaveh/CrackFlowFer2010/",recursive = TRUE)
file.copy(list.files("~/multibasinpaper/LOC_Ferrisburg_2015/",full.names = TRUE,recursive = TRUE),"/dev/shm/rojakaveh/CrackFlowFer2015/",recursive = TRUE)
####load updated functions
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/SWATmodel/R/readSWAT.R?root=ecohydrology")
save(readSWAT,file="readSWAT.R")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/setup_swatcal.R?root=ecohydrology")
source("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/EcoHydRology/R/swat_objective_function_rch.R?root=ecohydrology")

########load calib params
change_params=""
rm(change_params)
load(paste(path.package("EcoHydRology"), "data/change_params.rda", sep = "/"))
calib_range=c("1999-12-31","2021-12-31")
params_select=c(1,2,3,4,5,6,7,8,9,10,11,14,19,21,23,24,32,33)
calib_params=change_params[params_select,]
#calib_params$current[2]=0.5
calib_params$min[9]=0
calib_params$min[10]=0
calib_params$current[9]=2.5
calib_params$current[10]=2.5
calib_params$min[11]=0.01
calib_params$max[11]=1
calib_params$min[13]=40
calib_params$max[13]=95
calib_params$min[14]=0.3
calib_params$max[14]=3
#calib_params$min[15]=0.3
#calib_params$max[15]=3
calib_params$min[15]=0.3
calib_params$max[15]=3
calib_params$min[16]=0.3
calib_params$max[16]=3
calib_params$max[2]=2
calib_params$max[3]=600
calib_params$max[4]=0.3
calib_params$min[2]=0
calib_params$max[2]=1
calib_params$current[2]=0.5
calib_params[c(9,10),4]=0
calib_params[c(9,10),6]=2.5
calib_params[11,5]=1
calib_params[11,6]=.5
calib_params$min[1]=0


calib_params[1:7]
setup_swatcal(calib_params)
rch=3

# Test calibration
x=calib_params$current
swat_objective_function_rch=function (x, calib_range, calib_params, flowgage, rch,save_results=F)
{
  pacman::p_load(SWATmodel,dplyr,EcoHydRology,base,topmodel,utils)
  calib_params$current <- x
  tmpdir=as.character(as.integer((runif(1)+1)*10000))
  tmpdir=paste(c(format(Sys.time(), "%s"),tmpdir,Sys.getpid()),sep="",collapse="")
  print(tmpdir)
  dir.create(tmpdir)
  file.copy(list.files(),tmpdir)
  setwd(tmpdir)
  file.remove(list.files(pattern="output."))
  alter_files(calib_params)
  libarch = if (nzchar(base::version$arch)) paste("libs", base::version$arch, sep = "/") else "libs"
  swatbin <- "rswat2012.exe"
  junkout=system(shQuote(paste(path.package("SWATmodel"), libarch, swatbin, sep = "/")),intern = TRUE,ignore.stderr=TRUE)
  start_year = read.fortran(textConnection(readLines("file.cio")[9]), "f20")
  load("readSWAT.R")
  outdata = readSWAT("rch",".")
  test2 = subset(outdata, outdata$RCH == rch)
  test3 = merge(flowgage$flowdata, test2, all = F)
  NS = NSeff(test3$Qm3ps, test3$FLOW_OUTcms)
  print(NS)
  if(save_results){
    SWAToutput()
    for(diffname in grep("output|unixorig",grep(paste0(unique(calib_params[,1]),
                                                       collapse = "|"),list.files(),value=TRUE),
                         invert = TRUE,value=TRUE)){
      if(length(grep("calib",diffname))>0){next}
      junk1=readLines(diffname)
      junk2=readLines(paste0("../",diffname))
      difflocs=as.numeric(strsplit(paste0(as.data.frame(strsplit(ses(junk1,junk2),split = "c"))[2,],collapse = ","),split=",")[[1]])
      if(length(difflocs)>0){
        write(diffname,"calibdiffs.txt",append=TRUE)
        write(paste0(diffname,":",junk1[difflocs]," < ",junk2[difflocs]),"calibdiffs.txt",append=TRUE)
      }
    }
    calibdir=paste0("../calib",format(Sys.time(),format="%Y%m%d%H%M"))
    dir.create(calibdir)
    file.copy(list.files(),calibdir)
  }
  file.remove(list.files())
  setwd("../")
  file.remove(tmpdir)
  return(abs(NS - 1))
}
swat_objective_function_rch(x,calib_range,calib_params,flowgage,rch,save_results=F)

#################################
#Loop to have 10datasets of outdeoptim
for (calitts in 1:5){
  print(calitts)
  set.seed(calitts)
  cl <- parallel::makeCluster(16)
  outDEoptim<-DEoptim(swat_objective_function_rch,calib_params$min,calib_params$max,
                      DEoptim.control(cluster=cl,strategy = 6,NP = 16,itermax=300,parallelType = 1,
                                      packages = c("SWATmodel","dplyr","EcoHydRology","base","topmodel","utils"),parVar=c("%<%","NSeff","read.fortran","readSWAT","alter_files")),calib_range,calib_params,flowgage,rch)
  x=outDEoptim$optim$bestmem
  filename=format(Sys.time(), "%Y%m%d%H%M_outDEoptim.RData")
  assign(filename,outDEoptim)
  save(outDEoptim,file=format(Sys.time(), "~/multibasinpaper/%Y%m%d%H%M_outDEoptim.RData"))
}
