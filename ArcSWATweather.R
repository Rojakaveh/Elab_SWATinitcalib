source("https://raw.githubusercontent.com/Rojakaveh/FillMissWX/main/FillMissWX.R")#loading fillmisswx func
pacman::p_load(rnoaa,SWATmodel)
##-----getting usgs info--------------------------
flowgage_id="04282650" #Little Otter Creek at Ferrisburg, VT.
flowgage=get_usgs_gage(flowgage_id,begin_date = "2010-01-01",end_date= "2021-12-31")
flowgage$flowdata$Qmm=(flowgage$flowdata$flow)/(flowgage$area*1000)#mm
##getting ghcn weather data
WXData=FillMissWX(declat = flowgage$declat,declon = flowgage$declon,StnRadius =30,date_min=min(flowgage$flowdata$mdate),date_max=max(flowgage$flowdata$mdate),method = "IDW",minstns = 3,alfa = 2)
AllDays=data.frame(date=seq(min(flowgage$flowdata$mdate), by = "day", length.out = max(flowgage$flowdata$mdate)-min(flowgage$flowdata$mdate)))
WXData=merge(AllDays,WXData,all=T)
WXData$PRECIP=WXData$P
WXData$PRECIP[is.na(WXData$PRECIP)]=-99
WXData$TMX=WXData$MaxTemp
WXData$TMX[is.na(WXData$TMX)]=-99
WXData$TMN=WXData$MinTemp
WXData$TMN[is.na(WXData$TMN)]=-99
WXData$DATE=WXData$date
Roj0pcp=data.frame(matrix(0,ncol = 1,nrow = nrow(WXData)+1))
Roj0pcp[1,1]=paste0(year(WXData$DATE[1]),format.Date(WXData$DATE[1],"%d"),format.Date(WXData$DATE[1],"%m"))
Roj0pcp[2:nrow(Roj0pcp),]=sprintf("%.3f", round(WXData$PRECIP,3))
write.table(Roj0pcp,"Roj0pcp.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
pcpfork=matrix(0,nrow = 2,ncol = 1)
pcpfork[1,]=paste0("ID,NAME,LAT,LONG,ELEVATION")
pcpfork[2,]=paste0("1,","Roj0pcp,",format(round(flowgage$declat,3)),",",format(round(flowgage$declon,3)),",",format(round(flowgage$elev,3)))
write.table(pcpfork,"pcpfork.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
####
Roj0tmp=data.frame(matrix(0,ncol = 1,nrow = nrow(WXData)+1))
Roj0tmp[1,1]=paste0(year(WXData$DATE[1]),format.Date(WXData$DATE[1],"%d"),format.Date(WXData$DATE[1],"%m"))
Roj0tmp[2:nrow(Roj0tmp),]=paste0(sprintf("%.3f", round(WXData$TMX,3)),",",sprintf("%.3f", round(WXData$TMN,3)))
write.table(Roj0tmp,"Roj0tmp.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
tmpfork=matrix(0,nrow = 2,ncol = 1)
tmpfork[1,]=paste0("ID,NAME,LAT,LONG,ELEVATION")
tmpfork[2,]=paste0("1,","Roj0tmp,",format(round(flowgage$declat,3)),",",format(round(flowgage$declon,3)),",",format(round(flowgage$elev,3)))
write.table(tmpfork,"tmpfork.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

