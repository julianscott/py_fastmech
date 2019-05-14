
library(tidyverse)
## for three iterables
q_text <- "q1_8"
tmpwd <- paste0("E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\",q_text,"\\")
r1 <- read.csv(paste0(tmpwd,"case\\",q_text,"_penz_m_vs_s.csv"),header = FALSE)
head(r1)
str(r1)
r2 <- read.csv(paste0(tmpwd,"case\\",q_text,"_penz_m_rmse.csv"),header = FALSE)
colnames(r2) <- c("Cd_bs","Cd_dl","Cd_bc","rmse","meon")

Cdhead_bs <- paste0("Cd_bs",as.character(r2[1]))
Cdhead_bs <- sub("","Cd_bs",r2[[1]])
# Cdhead_bs
Cdhead_dl <- paste0("Cd_dl",as.character(r2[2]))
Cdhead_dl <- sub("","Cd_dl",r2[[2]])
# Cdhead_dl
Cdhead_bc <- paste0("Cd_bc",as.character(r2[3]))
Cdhead_bc <- sub("","Cd_bc",r2[[3]])
# Cdhead_bc
Cdhead <- paste0(Cdhead_bs,"/",Cdhead_dl,"/",Cdhead_bc)
colnames(r1) <- c("long","X","Y","ObsZ",Cdhead)
r3 <- r1 %>%
  gather(key = "Cd",value = "Z",-c(long,X,Y)) %>%
  mutate(long100m = long - long %% 100)

head(r3)
head(r2)
r2 <- arrange(r2,rmse) %>%
  mutate(Cdhead = paste0("Cd_bs",Cd_bs,"/","Cd_dl",Cd_dl,"/","Cd_bc",Cd_bc))
ggplot(data = r2) +
  geom_point(aes(x=Cdhead,y = rmse))
View(r2)
qlist <- c("Cd_bs0.003/Cd_dl0.008/Cd_bc0.0166666666666667")

ggplot() + 
  geom_line(data = filter(r3,Cd == "ObsZ"),aes(x = long,y = Z),color = "grey",size = 2) +
  geom_line(data = filter(r3,Cd %in% qlist),aes(x = long,y = Z,color = Cd),size = 1) +
  geom_point(data = filter(r3,Cd == "ObsZ"),aes(x = long,y = Z),fill = "grey",color = "black") +
  geom_point(data = filter(r3,Cd %in% qlist),aes(x = long,y = Z,fill = Cd,color = Cd)) + 
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
# facet_wrap(long100m ~. ,ncol = 1,scales = "free")
