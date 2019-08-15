# see https://www.gastonsanchez.com/r4strings/input-and-output.html

# If “.” matches any character, how do you match a literal “.”? You need to use an “escape” to 
# tell the regular expression you want to match it exactly, not use its special behaviour. 
# Like strings, regexps use the backslash, \, to escape special behaviour. So to match an ., 
# you need the regexp \.. Unfortunately this creates a problem. We use strings to represent 
# regular expressions, and \ is also used as an escape symbol in strings. So to create the 
# regular expression \. we need the string "\\.".

# Discharge
q = 400
### downstream boundary condition WSE 
H_DS =   78.28746

basedir = "E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\"
solverdir = "C:\\Users\\jascott\\iRIC\\solvers\\fastmech\\"
study_area_coords = "E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\Demo_Directory\\Multiple_Discharge_Code\\q3_8_to_8_0\\smrf_DEM_v24_points_penz_m.txt"

# drag coefficient for class 0 base sand
cdmin_Cd0 = 0.005
cdmax_Cd0 = 0.009
cdn_Cd0 = 2

# drag coefficient for class 1 dogleg
cdmin_Cd1 = 0.05
cdmax_Cd1 = 0.09
cdn_Cd1 = 2

# drag coefficient for class 2 bld/cob
cdmin_Cd2 = 0.02
cdmax_Cd2 = 0.035
cdn_Cd2 = 2

# drag coefficient for class 3 veg
cdmin_Cd3 = 0.005
cdmax_Cd3 = 0.009
cdn_Cd3= 3

# drag coefficient for class 4
cdmin_Cd4 = 0.02
cdmax_Cd4 = 0.02
cdn_Cd4= 1

# 1D step back water intial condition = 2
iniType = 2
# 1D initial condition drag coef (controls to intial WSE slope)
OneDCD = 0.01

# LEV switch on or off
LEV_Type = 0
LEV_Constant = 0.06
LEVStart_Iter = 500
LEVEnd_Iter = 1000
StartLEV = 0.06
EndLEV = 0.006

qdir = if_else(as.integer(q) == q, paste0("q",q,"_0"),
               paste0("q",sub(pattern = "\\.",replacement = "_",x = as.character(q))))
qdir2 = if_else(as.integer(q) == q, paste0(q,"_0"),
                paste0(sub(pattern = "\\.",replacement = "_",x = as.character(q))))


write.table("", paste0(basedir,qdir,"\\",qdir,"cms_config.ini"),quote = F,row.names = F,col.names = F)

inifile <- paste0(basedir,qdir,"\\",qdir,"cms_config.ini")

cat("[Params]","\n",file = inifile)
cat("meas_WSE_File = ",basedir,qdir,"\\","ReadObsWSE_XYZ_",qdir2,".csv","\n",sep="",file = inifile,append = T) 
cat("study_area_coords =", study_area_coords,"\n","\n",sep="",file = inifile,append = T)

cat("# drag coefficient for class 0 base sand","\n",file = inifile,append = T)
cat("cdmin_Cd0 = ", cdmin_Cd0,"\n",file = inifile,append = T)
cat("cdmax_Cd0 = ", cdmax_Cd0,"\n",file = inifile,append = T)
cat("cdn_Cd0 = ", cdn_Cd0,"\n","\n",file = inifile,append = T)

cat("# drag coefficient for class 1 dogleg","\n",file = inifile,append = T)
cat("cdmin_Cd1 = ",cdmin_Cd1,"\n",file = inifile,append = T)
cat("cdmax_Cd1 = ",cdmax_Cd1,"\n",file = inifile,append = T)
cat("cdn_Cd1 = ",cdn_Cd1,"\n","\n",file = inifile,append = T)

cat("# drag coefficient for class 2 bld/cob")
cat("cdmin_Cd2 = ",cdmin_Cd2,"\n",file = inifile,append = T)
cat("cdmax_Cd2 = ",cdmax_Cd2,"\n",file = inifile,append = T)
cat("cdn_Cd2 = ",cdn_Cd2,"\n","\n",file = inifile,append = T)

cat("# drag coefficient for class 3 veg","\n",file = inifile,append = T)
cat("cdmin_Cd3 = ",cdmin_Cd3,"\n",file = inifile,append = T)
cat("cdmax_Cd3 = ",cdmax_Cd3,"\n",file = inifile,append = T)
cat("cdn_Cd3 = ",cdn_Cd3,"\n","\n",file = inifile,append = T)

cat("# drag coefficient for class 4","\n",file = inifile,append = T)
cat("cdmin_Cd4 = ",cdmin_Cd4,"\n",file = inifile,append = T)
cat("cdmax_Cd4 = ",cdmax_Cd4,"\n",file = inifile,append = T)
cat("cdn_Cd4 = ",cdn_Cd4,"\n","\n",file = inifile,append = T)

cat("Q = ",q,"\n",file = inifile,append = T)
cat("H_DS = ",H_DS,"\n",file = inifile,append = T)
cat("LEV_type = ", LEV_Type,"\n",file = inifile,append = T)

cat("LEV_Constant = ", LEV_Constant,"\n",file = inifile,append = T)
cat("LEVStart_Iter = ", LEVStart_Iter,"\n",file = inifile,append = T)
cat("LEVEnd_Iter = ", LEVEnd_Iter,"\n",file = inifile,append = T)
cat("StartLEV = ", StartLEV,"\n",file = inifile,append = T)
cat("EndLEV = ", EndLEV,"\n",file = inifile,append = T)

cat("iniType =", iniType,"\n",file = inifile,append = T)
cat("OneDCD =", OneDCD,"\n",file = inifile,append = T)
cat("working_dir =",basedir,qdir,"\\case","\n","\n",sep="",file = inifile,append = T)

cat("solver_path = ;",solverdir,"\n",sep="",file = inifile,append = T)

cat("base_file = ",basedir,qdir,"\\case\\Case1.cgn","\n",sep="",file = inifile,append = T)
cat("rmse_file = q570_0_penz_m_rmse.csv","\n",file = inifile,append = T)
cat("meas_vs_sim_file = q570_0_penz_m_vs_s.csv","\n",file = inifile,append = T)

