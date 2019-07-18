# load rating curves
# As of 03/29/19
library(tidyverse)
library(segmented)

###############################################################
### Build or load rating curves

uprc <- read.csv("E:\\_DoD\\_Camp_Pendleton_Survey\\Hydrology\\Pressure_Transducer_data\\PT_data_Time_Flow_Alignment_work\\031219_fullrecord\\Upper_PT\\Upper_PT_RatingCurve_Table.csv")
mdrc <- read.csv("E:\\_DoD\\_Camp_Pendleton_Survey\\Hydrology\\Pressure_Transducer_data\\PT_data_Time_Flow_Alignment_work\\031219_fullrecord\\Middle_PT\\Middle_PT_RatingCurve_Table.csv")
dnrc <- read.csv("E:\\_DoD\\_Camp_Pendleton_Survey\\Hydrology\\Pressure_Transducer_data\\PT_data_Time_Flow_Alignment_work\\031219_fullrecord\\Down_PT\\rating_curve_table_DownStream_PT.csv")

uprc <- arrange(uprc,q) %>%
  filter( X != 3) # removed because not monotonic and removal is inconsequential.
mdrc <- arrange(mdrc,q)
dnrc <- arrange(dnrc,q)

# test to make sure series is monotonic
dnrc %>%
  mutate(hlag = lag(h)) %>%
  summarize(min_delta = min(h-hlag,na.rm = T))


###############################################################
### Function for making predictions. Arguments are a rating table and a discharge value
### Tweak to use predict() for a nls model, segmented model, etc

# get stage for a given rating table
rc_approx<- function(rc_tbl,q){
  h = approx(x=rc_tbl$q,y=rc_tbl$h,xout = q)$y
  return(h)
}

###############################################################
### Define some discharges to get stages from
Q_calib <- seq(2.7,2.9,0.1)
Q_calib <- 33
Q_calib
###############################################################
### Use function to retrieve stage values 
up_h <- rc_approx(uprc,Q_calib)
mid_h <- rc_approx(mdrc,Q_calib)
dwn_h <- rc_approx(dnrc,Q_calib)
dwnext_h <- rc_approx(dwnext_rc,Q_calib)

###############################################################
### Build csv files for fastmech/python

### build ObsWSE csv files (one per discharge, with PT coords + WSE) 
uppt = c(1033718.763,3709866.074)
midpt = c(1033545.146,3710054.03)
dwnpt = c(1033528.961,3710106.519)
ptdf = data.frame("X" = numeric(0),"Y" = numeric(0))
ptdf[1,] = uppt
ptdf[2,] = midpt
ptdf[3,] = dwnpt
ptdf$pt = c("up","mid","dwn")
ptdf$long = c(1.5,241.72,294.54) # these longitudinal posiotions are unique to the SMR
ptdf <- ptdf %>%
  dplyr::select(pt,long,X,Y) 

# directory
q_text <- "q8_0"
dir1 <- paste0("E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\",q_text,"\\")
# list.files(dir1)
# i=1
# For Roughness Varying Python Tool
for(i in 1:length(Q_calib)){
  dfi = ptdf 
  q = Q_calib[i]
  dfi$ObsWSE = sapply(list(uprc,mdrc,dnrc),function(x) rc_approx(rc_tbl=x,q=q))
  q = ifelse(q %% 1 == 0,paste0(q,"_0"),sub("\\.","_",as.character(q)))
  dfgui = dfi %>%
    dplyr::select(X,Y,ObsWSE)
  colnames(dfgui)[which(colnames(dfgui) == "ObsWSE")] <- q
  write.table(dfgui,paste0(dir1,"ReadObsWSE_XYZ_",q,"_gui.csv"),sep = ",",qmethod = "double",col.names = T,row.names = F)
  dfpy= dfi %>%
    dplyr::select(pt,long,X,Y,ObsWSE)
  colnames(dfpy)[which(colnames(dfpy) == "ObsWSE")] <- q
  write.table(dfpy,paste0(dir1,"ReadObsWSE_XYZ_",q,".csv"),sep = ",",qmethod = "double",col.names = T,row.names = F)
}

# # For Discharge Batch Process Python Code
# for(i in 1:length(Q_calib)){
#   dfi = ptdf 
#   q = Q_calib[i]
#   dfi$ObsWSE = sapply(list(uprc,mdrc,dnrc),function(x) rc_approx(rc_tbl=x,q=q))
#   q = ifelse(q %% 1 == 0,paste0(q,"_0"),sub("\\.","_",as.character(q)))
#   dfgui = dfi %>%
#     dplyr::select(X,Y,ObsWSE)
#   colnames(dfgui)[which(colnames(dfgui) == "ObsWSE")] <- q
#   write.table(dfgui,paste0(dir1,"ReadObsWSE_XYZ_",q,".csv"),sep = ",",qmethod = "double",col.names = T,row.names = F)
# }
# 
# # For Discharge Batch Process v2 Python Code
# param_df <- data.frame(Q = character(length(Q_calib))) %>%
#   mutate(Q = Q_calib,
#          H_DS = rc_approx(dwnext_rc,Q))
# write.table(param_df,paste0(dir1,"QHDS_read.csv"),sep = ",",col.names = F,row.names = F)
