# Code for running fastmech models with varying roughness values, where three classes of roughness nodes are set to vary.
# The code copies the base file (Case1.cgn) that is PRECONFIGURED by the user, with roughness mapped to three classes of 
# roughness. Classes are 0, 1, and 5. Code will iteritevly re-classify these roughness patches for defined ranges of Cd.
# Tool is developed and annotated by Julian Scott, but is strongly dependent on code written by Rich McDonald, especially the VTK, h5py 
# utilization, and fastmech execution.

# 04/10/19

from __future__ import print_function  # Only Python 2.x
import numpy as np
from shutil import copyfile
import h5py
import vtk
import subprocess
import os
from itertools import count
import configparser
import re
import pandas as pd
import datetime

# Remember: any observed WSE points that do not fall onto the iric mesh will cause an error. They must not be included.

def getCellValue(vtkSGrid2D, newPoint2D, cellID, valarray):
    pcoords = [0.0,0.0,0.0]
    weights = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    clspoint = [0.,0.,0.]
    tmpid = vtk.mutable(0)
    vtkid2 = vtk.mutable(0)
    vtkcell2D = vtk.vtkQuad()
    vtkcell2D = vtkSGrid2D.GetCell(cellID)
    tmpres = vtkcell2D.EvaluatePosition(newPoint2D,clspoint,tmpid,pcoords,vtkid2, weights )
#    print(newPoint2D, clspoint, tmpid, pcoords, vtkid2, weights)
    idlist1 = vtk.vtkIdList()
    numpts = vtkcell2D.GetNumberOfPoints()
    idlist1 = vtkcell2D.GetPointIds()
    tmpVal = 0.0
    for x in range(0,numpts):
        tmpVal = tmpVal + weights[x]*valarray.GetTuple(idlist1.GetId(x))[0]
    return tmpVal

def isCellWet(vtkSGrid2D, newPoint2D, cellID, IBC_2D):
    pcoords = [0.0,0.0,0.0]
    weights = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    clspoint = [0.,0.,0.]
    tmpid = vtk.mutable(0)
    vtkid2 = vtk.mutable(0)
    vtkcell2D = vtk.vtkQuad()
    vtkcell2D = vtkSGrid2D.GetCell(cellID)
    tmpres = vtkcell2D.EvaluatePosition(newPoint2D,clspoint,tmpid,pcoords,vtkid2, weights )
    idlist1 = vtk.vtkIdList()
    numpts = vtkcell2D.GetNumberOfPoints()
    idlist1 = vtkcell2D.GetPointIds()
    tmpIBC = 0.0
    for x in range(0,numpts):
        tmpIBC = tmpIBC + weights[x]*abs(IBC_2D.GetTuple(idlist1.GetId(x))[0])
        # print(tmpIBC,abs(IBC_2D.GetTuple(idlist1.GetId(x))[0]))
    if tmpIBC >= .9999999:
        return True
    else:
        return False

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

def gen_filenames(prefix, suffix, places=3):
    """Generate sequential filenames with the format <prefix><index><suffix>

       The index field is padded with leading zeroes to the specified number of places

       http://stackoverflow.com/questions/5068461/how-do-you-increment-file-name-in-python
    """
    pattern = "{}{{:0{}d}}{}".format(prefix, places, suffix)
    for i in count(1):
        yield pattern.format(i)

def fastmech_change_cd(hdf_file, newCd):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttCD/Value']
    dset = group[u' data']
    # print dset[0]
    dset[0] = newCd
    # print dset[0]
    file.close()

def fastmech_BCs(hdf_file, Q, H_DS, iniType, OneDCD, Cddistribution):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttQ/Value']
    dset = group[u' data']
    dset[0] = Q
    group2 = file['/iRIC/CalculationConditions/FM_HydAttCDType/Value']
    dset2 = group2[u' data']
    dset2[0] = Cddistribution
    group3 = file['/iRIC/CalculationConditions/FM_HydAttWS/Value']
    dset3 = group3[u' data']
    dset3[0] = H_DS
    group4 = file['/iRIC/CalculationConditions/FM_HydAttWSType/Value']
    dset4 = group4[u' data']
    dset4[0] = iniType
    # group5 = file['/iRIC/CalculationConditions/FM_HydAttWS1DStage/Value']
    # dset5 = group5[u' data']
    # dset5[0] = OneDStage
    # group6 = file['/iRIC/CalculationConditions/FM_HydAttWS1DDisch/Value']
    # dset6 = group6[u' data']
    # dset6[0] = OneDQ
    group7 = file['/iRIC/CalculationConditions/FM_HydAttWS1DCD/Value']
    dset7 = group7[u' data']
    dset7[0] = OneDCD
    file.close()


#hdf_file = hdf5_file_name
#newCd_0 = Cd
#newCd_1 = 0.0032

# for changing constant roughness for nodes set to 0 (i.e. vary 1 class of roughness nodes)    
def fastmech_change_var_cd(hdf_file, newCd):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/iRICZone/GridConditions/roughness/Value']
    dset = group[u' data']
    for index, val in enumerate(dset):
        if val == 0.0:
            dset[index] = newCd
#         else:
#            dset[index] = newCd_1 #keep values in original project, change only values with 0
    # print dset[0]
    # print dset[0]
    file.close()
    
#hdf_file = hdf5_file_name
#newCd_Cd0 = Cd_Cd0
#newCd_dl = Cd_Cd1
#newCd_Cd5 = Cd_Cd5

# For changing roughness for nodes set to 0,1,and 5 (i.e. vary 3 class of roughness nodes)  
def fastmech_change_var_cd2(hdf_file, newCd_Cd0,newCd_Cd1,newCd_Cd5):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/iRICZone/GridConditions/roughness/Value']
    dset = group[u' data']
    for value in list([0.0,1.0,5.0]):
        if value in dset:
            print("")
        else:
            return(print("There are not the correct three classes of roughness in the base file"))
    for index, val in enumerate(dset):
        if val == 0.0:
            dset[index] = newCd_Cd0 # base cd class
        else:
            if val == 1.0:
                dset[index] = newCd_Cd1 # dogleg cd class
            else:
                if val == 5.0:
                    dset[index] = newCd_Cd5 # boulder/cobble cd class
    file.close()

def check_Cd_classes(hdf_file,Cd0_class,Cd1_class,Cd5_class):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/iRICZone/GridConditions/roughness/Value']
    dset = group[u' data']
    for value in list([Cd0_class,Cd1_class,Cd5_class]):
        if value in dset:
            print("Roughness class " + str(value) + " confirmed to be present")
        else:
            return(print("There are not the correct three classes of roughness in the base file"))
    file.close()


#hdf_file = hdf5_file_name
#newWSE_0 = H_DSmin
#newWSE_1 = 0.0032
def fastmech_change_var_H_DS(hdf_file, newWSE_0):
    # hdf5_file_name = r'F:\Kootenai Project\USACE\Braided\Case11_tmp.cgn'
    # r+ adds read/write permisions to file
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttWS/Value']
    dset = group[u' data']
    for index, val in enumerate(dset):
        if val == 0.0:
            dset[index] = newWSE_0
    file.close()

# Update other model parameters from variables defined in the config file.
def fastmech_params(hdf_file, Itermax):
    file = h5py.File(hdf_file, 'r+')
#    group = file['/iRIC/CalculationConditions/FM_EndLEV/Value']
#    dset = group[u' data']
#    dset[0] = endLEV
    group = file['/iRIC/CalculationConditions/FM_SolAttItm/Value']
    dset = group[u' data']
    dset[0] = Itermax
    file.close()

def add_fastmech_solver_to_path(solverPath):
    os.environ['PATH'] += solverPath

def create_vtk_structured_grid(sgrid, hdf5_file_name, xoffset, yoffset):
    # type: (object) -> object
    file = h5py.File(hdf5_file_name, 'r')
    xcoord_grp = file['/iRIC/iRICZone/GridCoordinates/CoordinateX']
    # print(xcoord_grp.keys())
    ycoord_grp = file['/iRIC/iRICZone/GridCoordinates/CoordinateY']
    # print(ycoord_grp.keys())
    wse_grp = file['iRIC/iRICZone/FlowSolution1/WaterSurfaceElevation']
    # print(wse_grp.keys())
    topo_grp = file['iRIC/iRICZone/FlowSolution1/Elevation']
    # print(topo_grp.keys())
#    ibc_grp = file['iRIC/iRICZone/FlowSolution1/IBC']
#    velx_grp = file['iRIC/iRICZone/FlowSolution1/VelocityX']
#    vely_grp = file['iRIC/iRICZone/FlowSolution1/VelocityY']
    depth_grp = file['iRIC/iRICZone/FlowSolution1/Depth']

    xcoord_data = xcoord_grp[u' data']
    ycoord_data = ycoord_grp[u' data']
    wse_data = wse_grp[u' data']
    topo_data = topo_grp[u' data']
#    ibc_data = ibc_grp[u' data']
#    velx_data = velx_grp[u' data']
#    vely_data = vely_grp[u' data']
    depth_data = depth_grp[u' data']
    # SGrid = vtk.vtkStructuredGrid()
    ny, nx, = xcoord_data.shape
    # print(ny, nx)
    sgrid.SetDimensions(nx, ny, 1)
    points = vtk.vtkPoints()
    wseVal = vtk.vtkFloatArray()
    wseVal.SetNumberOfComponents(1)
#    ibcVal = vtk.vtkIntArray()
#    ibcVal.SetNumberOfComponents(1)
#    velVal = vtk.vtkFloatArray()
#    velVal.SetNumberOfComponents(1)
    depthVal = vtk.vtkFloatArray()
    depthVal.SetNumberOfComponents(1)
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(xcoord_data[j, i] - xoffset, ycoord_data[j, i] - yoffset, 0.0)
            wseVal.InsertNextValue(wse_data[j, i])
#            ibcVal.InsertNextValue(ibc_data[j, i])
#            velVal.InsertNextValue(np.sqrt(np.power(velx_data[j, i],2) + np.power(vely_data[j,i],2)))
            depthVal.InsertNextValue(depth_data[j, i])
        sgrid.SetPoints(points)

        sgrid.GetPointData().AddArray(wseVal)
#        sgrid.GetPointData().AddArray(ibcVal)
#        sgrid.GetPointData().AddArray(velVal)
        sgrid.GetPointData().AddArray(depthVal)
    wseVal.SetName("WSE")
#    ibcVal.SetName("IBC")
#    velVal.SetName("Velocity")
    depthVal.SetName("Depth")
    
# Same as above, but also calcs velocity for final calculation of mean depth and mean velocity
def create_vtk_structured_grid_post(sgrid, hdf5_file_name, xoffset, yoffset):
    # type: (object) -> object
    file = h5py.File(hdf5_file_name, 'r')
    xcoord_grp = file['/iRIC/iRICZone/GridCoordinates/CoordinateX']
    # print(xcoord_grp.keys())
    ycoord_grp = file['/iRIC/iRICZone/GridCoordinates/CoordinateY']
    # print(ycoord_grp.keys())
    wse_grp = file['iRIC/iRICZone/FlowSolution1/WaterSurfaceElevation']
    # print(wse_grp.keys())
    topo_grp = file['iRIC/iRICZone/FlowSolution1/Elevation']
    # print(topo_grp.keys())
#    ibc_grp = file['iRIC/iRICZone/FlowSolution1/IBC']
    velx_grp = file['iRIC/iRICZone/FlowSolution1/VelocityX']
    vely_grp = file['iRIC/iRICZone/FlowSolution1/VelocityY']
    depth_grp = file['iRIC/iRICZone/FlowSolution1/Depth']

    xcoord_data = xcoord_grp[u' data']
    ycoord_data = ycoord_grp[u' data']
    wse_data = wse_grp[u' data']
    topo_data = topo_grp[u' data']
#    ibc_data = ibc_grp[u' data']
    velx_data = velx_grp[u' data']
    vely_data = vely_grp[u' data']
    depth_data = depth_grp[u' data']
    # SGrid = vtk.vtkStructuredGrid()
    ny, nx, = xcoord_data.shape
    # print(ny, nx)
    sgrid.SetDimensions(nx, ny, 1)
    points = vtk.vtkPoints()
    wseVal = vtk.vtkFloatArray()
    wseVal.SetNumberOfComponents(1)
#    ibcVal = vtk.vtkIntArray()
#    ibcVal.SetNumberOfComponents(1)
    velVal = vtk.vtkFloatArray()
    velVal.SetNumberOfComponents(1)
    depthVal = vtk.vtkFloatArray()
    depthVal.SetNumberOfComponents(1)
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(xcoord_data[j, i] - xoffset, ycoord_data[j, i] - yoffset, 0.0)
            wseVal.InsertNextValue(wse_data[j, i])
#            ibcVal.InsertNextValue(ibc_data[j, i])
            velVal.InsertNextValue(np.sqrt(np.power(velx_data[j, i],2) + np.power(vely_data[j,i],2)))
            depthVal.InsertNextValue(depth_data[j, i])
        sgrid.SetPoints(points)

        sgrid.GetPointData().AddArray(wseVal)
#        sgrid.GetPointData().AddArray(ibcVal)
        sgrid.GetPointData().AddArray(velVal)
        sgrid.GetPointData().AddArray(depthVal)
    wseVal.SetName("WSE")
#    ibcVal.SetName("IBC")
    velVal.SetName("Velocity")
    depthVal.SetName("Depth")


# Set the config file 
setfile =  "E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\q8_0\\q8_0cms_config.ini"
# Set up configuration file parser
config = configparser.ConfigParser()
# Read in the configuration file
config.read(setfile)

# read in coordinates of study area grid cell centers
# csv format: pt number,X,Y,Z
#study_area_coords = "E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\smrf_DEM_v24_points_penz_m.txt"
#study_area_tbl = np.genfromtxt(study_area_coords, delimiter=',', skip_header=0)
study_area_tbl = np.genfromtxt(config.get('Params','study_area_coords'), delimiter=',', skip_header=0)

# Get the directory location of the Obs WSE csv file (must be formatted correctly)
meas_wse = pd.read_csv(config.get('Params','meas_WSE_File'))
# Convert to a numpy array
meas_wse = np.array(meas_wse)
# Number of Obs WSEs
nummeas = meas_wse.shape[0]

# Temporarily add fastmech solver to environment paths
add_fastmech_solver_to_path(config.get('Params','lib_path'))
add_fastmech_solver_to_path(config.get('Params','solver_path'))

# Set working directory from parameter in config file
os.chdir(config.get('Params','working_dir'))

# Set naming convention for the produced fastmech CGNS file 
# This option, with other tweaks to the code below, can be used to produce a new file for each iteration (However, for many runs = many gigabytes)
#g = gen_filenames("FM_Calib_Flow_", ".cgns") 
# This configuraton of the code produces a CGNS file that overwrites itself each iteration. 
# Then, just the best model is written at the conclusion of all iterations.
g = "FM_Calib_Flow2.cgns" 

# Configured for 3 roughness classes
# Cd0 = base nodes

cdmin_Cd0 = eval(config.get('Params','cdmin_Cd0'))
cdmax_Cd0 = eval(config.get('Params','cdmax_Cd0'))
cdn_Cd0 = eval(config.get('Params','cdn_Cd0'))
# Cd1 = dogleg reach nodes
cdmin_Cd1 = eval(config.get('Params','cdmin_Cd1'))
cdmax_Cd1 = eval(config.get('Params','cdmax_Cd1'))
cdn_Cd1 = eval(config.get('Params','cdn_Cd1'))
# Cd5 = boulder/cobble nodes
cdmin_Cd5 = eval(config.get('Params','cdmin_Cd5'))
cdmax_Cd5 = eval(config.get('Params','cdmax_Cd5'))
cdn_Cd5 = eval(config.get('Params','cdn_Cd5'))

# This is from Rich's original code; Leave set to 0.
#xoffset = config.getfloat('Params','xoffset')
#yoffset = config.getfloat('Params','yoffset')
xoffset = 0
yoffset = 0
# Discharge
Q = config.getfloat('Params','Q')
# Down stream WSE boundary condition
H_DS = config.getfloat('Params','H_DS')
#eval(config.get('Params','H_DS'))
# Roughness distribution type; val = 1 for roughness variable by node 
Cddistribution = 1 
# Initial condition type 
iniType = config.getfloat('Params','iniType')
# Drag coef for 1D Step backwater initial condition model
OneDCD = config.getfloat('Params','OneDCD')

# Solution iterations 
sol_iterations = 1500

# Create an array for each roughness values/class
a0 = np.linspace(cdmin_Cd0,cdmax_Cd0,int(cdn_Cd0))
a1 = np.linspace(cdmin_Cd1,cdmax_Cd1,int(cdn_Cd1))
a5 = np.linspace(cdmin_Cd5,cdmax_Cd5,int(cdn_Cd5))
# Set the number of roughness classes
n_cd_classes = 3

# Create an array of all possible combinations of roughness values/class
# by gridding arrays and reshaping into n_cd_classes columns
numcds = np.array(np.meshgrid(a0,a1,a5)).T.reshape(-1,n_cd_classes)
print("number of iterations: " + str(numcds.shape[0]))

# Create containers for rmse, mean error on discharge, and measured vs simulated output (1 addition per containter per iteration)
rmse_data = np.zeros(numcds.shape[0])
meod_val = np.zeros(numcds.shape[0])
meas_and_sim_wse = np.zeros(shape=(nummeas, numcds.shape[0]+1))
# Create containers for capturing current iteration cd values for inclusion in output
cd_val_Cd0 = np.zeros(numcds.shape[0])
cd_val_Cd1 = np.zeros(numcds.shape[0])
cd_val_Cd5 = np.zeros(numcds.shape[0])

# Check to make sure that the Cd classes (0,1, and 5) are present in the initial Case1.cgn file.
# If only one or two Cd classes are desired, omit all respective variables from the above code 
# and the code should be functional still, with minor tweaks I expect.
Cd0_class = 0 # base/sand/pebble
Cd1_class = 1 # dogleg
Cd5_class = 5 # boulder cobble

# Set a file name to the hdf5 file (this is the fastmech cgns file produced when the base_file (Case1.cgn) is copied)
hdf5_file_name = g
# Copy the base_file (Case1.cgn). 
# IMPORTANT: The resulting file inherits all grid values (elevation, roughness, etc) and Calculation Conditions. 
#            If a given Calculation Condition or grid node value is not modified by this code, it will therefore be
#            equal to whatever was inherited from the base_file.
copyfile(config.get('Params','base_file'),hdf5_file_name)
check_Cd_classes(hdf5_file_name,Cd0_class,Cd1_class,Cd5_class)

# Check how many WSE elevations are present
print(str(nummeas) + " WSE observations ready for calibration")
# Check statement to confirm important batchmode settings
print("Begin fastmech model runs for Q = ",str(Q),", iterating over", str(numcds.shape[0]) ,"combinations of drag coefficients")

# get start time for timing batch process
start = datetime.datetime.now()
cdind = 0
# Begin batch process`
# For each value(cdind) in the range 0:numcds:
for cdind in range(numcds.shape[0]):
    # Assign values to each Cd class, for the given cdind
    Cd_Cd0 = numcds[cdind,0]
    Cd_Cd1 = numcds[cdind,1]
    Cd_Cd5 = numcds[cdind,2]
    # For the given fastmech run, populate containers to track Cd value by class
    cd_val_Cd0[cdind] = Cd_Cd0
    cd_val_Cd1[cdind] = Cd_Cd1
    cd_val_Cd5[cdind] = Cd_Cd5

    # Use this option to create a new CGNS file for each iteration (many runs = many gigabytes)
    # hdf5_file_name = next(g)
    # Set CGNS name (will be overwritten)
    hdf5_file_name = g
    # Copy the base_file (Case1.cgn). 
    # IMPORTANT: The resulting file inherits all grid values (elevation, roughness, etc) and Calculation Conditions. 
    #            If a given Calculation Condition or grid node value is not modified by this code, it will therefore be
    #            equal to whatever was inherited from the base_file.
    copyfile(config.get('Params','base_file'),hdf5_file_name)
    # Update boundary condition values
    fastmech_BCs(hdf5_file_name, Q, H_DS, iniType, OneDCD, Cddistribution)
    # Update solution iterations. 
    fastmech_params(hdf5_file_name, sol_iterations)
    # Reclassify grid nodes for Cd roughness classes
    fastmech_change_var_cd2(hdf5_file_name, newCd_Cd0 = Cd_Cd0,newCd_Cd1 = Cd_Cd1,newCd_Cd5 = Cd_Cd5) 
    # Run fastmech solver using the configured hdf5 file
    for path in execute(["Fastmech.exe", hdf5_file_name]):
        try:
            meod = re.search('Discharge:(.+?)\n', path).group(1)  # For normal run
#            meod = re.search('(.+?)\n', path).group(1) # for LEV variable run
        except AttributeError:   
            meod = '-9999' # apply your error handling
#    Construct grid with solution values and create vtk modules or tools for extracting data by point
    SGrid = vtk.vtkStructuredGrid()
    try:
        create_vtk_structured_grid(SGrid, hdf5_file_name, xoffset, yoffset) 
    except KeyError:
#        When fastmech solution is not found, populate result boxes with -9999
        # Create container for simulated WSEs at the coordinates of the Obs WSEs    
        print('No Fastmech solution for run: ' + str(cdind))
        simwse = np.zeros(meas_wse.shape[0])
        # Create container for Obs WSEs
        measwse = np.zeros(meas_wse.shape[0])
        for counter, line in enumerate(meas_wse):
            point2D = [line[0+2]-xoffset, line[1+2]-yoffset, 0.0]
            pt1 = [line[0+2]-xoffset, line[1+2]-yoffset, 10.0]
            pt2 = [line[0+2]-xoffset, line[1+2]-yoffset, -10]
            # populate the simulated WSE container at position 'counter'
            simwse[counter] = -9999
            # populate the Observed WSE container at position 'counter'
            measwse[counter] = line[2+2]
             # For the given fastmech run, populate measured vs simulated WSE container with the simulated WSE data (array[nrow = 1,ncols = # ObsWSEs])
        meas_and_sim_wse[:,cdind] = simwse
        # For the given fastmech run, populate rmse container with results of rmse calculation (n = 1)
        rmse_data[cdind] = rmse(simwse, measwse)

        # For the given fastmech run, write the mean error on discharge
        meod_val[cdind] = -9999
        # For the given fastmech run, combine (sort of equivalent to R cbind) Cd class values, rmse, and mean error on discharge. 
        trmse = np.column_stack((cd_val_Cd0.flatten(),cd_val_Cd1.flatten(),cd_val_Cd5.flatten(),rmse_data.flatten(),meod_val.flatten()))
    else:
#        Else, when there is a Fastmech solution, create VTK grid and get results
        # The final string in the solution iteration output contains the mean error on discharge.
        meod =  "".join(meod.split())
        meod = float(meod)
        # Print Cd class values and mean error on discharge for the current model
        print("Cd_base:",Cd_Cd0)
        print("Cd_dogleg:",Cd_Cd1)
        print("Cd_boulder/cobble:",Cd_Cd5)
        print("Mean Error On Discharge:",meod)
        
        # Using the VTK grid created in line 462        
        cellLocator2D = vtk.vtkCellLocator()
        cellLocator2D.SetDataSet(SGrid)
        cellLocator2D.BuildLocator()
        WSE_2D = SGrid.GetPointData().GetScalars('WSE')
    #    IBC_2D = SGrid.GetPointData().GetScalars('IBC')
        Depth_2D = SGrid.GetPointData().GetScalars('Depth')
    #    Velocity_2D = SGrid.GetPointData().GetScalars('Velocity')
        
        # Create container for simulated WSEs at the coordinates of the Obs WSEs    
        simwse = np.zeros(meas_wse.shape[0])
        # Create container for Obs WSEs
        measwse = np.zeros(meas_wse.shape[0])
    #    counter = 0
    #    line = meas_wse[counter]
        # For each ObWSE point, extract the simulated WSE 
        for counter, line in enumerate(meas_wse):
            point2D = [line[0+2]-xoffset, line[1+2]-yoffset, 0.0]
            pt1 = [line[0+2]-xoffset, line[1+2]-yoffset, 10.0]
            pt2 = [line[0+2]-xoffset, line[1+2]-yoffset, -10]
            idlist1 = vtk.vtkIdList()
            cellLocator2D.FindCellsAlongLine(pt1,pt2,0.0, idlist1)
            cellid = idlist1.GetId(0)
            # Get WSE cell value
            tmpwse = getCellValue(SGrid, point2D, cellid, WSE_2D)
            # populate the simulated WSE container at position 'counter'
            simwse[counter] = tmpwse
            # populate the Observed WSE container at position 'counter'
            measwse[counter] = line[2+2]
            
        # For the given fastmech run, populate measured vs simulated WSE container with the simulated WSE data (array[nrow = 1,ncols = # ObsWSEs])
        meas_and_sim_wse[:,cdind] = simwse
        # For the given fastmech run, populate rmse container with results of rmse calculation (n = 1)
        rmse_data[cdind] = rmse(simwse, measwse)
        # Print RMSE
        print("rmse on Obs WSEs:",rmse_data[cdind])
#        # For the given fastmech run, populate containers to track Cd value by class
#        cd_val_Cd0[cdind] = Cd_Cd0
#        cd_val_Cd1[cdind] = Cd_Cd1
#        cd_val_Cd5[cdind] = Cd_Cd5
        # For the given fastmech run, write the mean error on discharge
        meod_val[cdind] = meod
        # For the given fastmech run, combine (sort of equivalent to R cbind) Cd class values, rmse, and mean error on discharge. 
        trmse = np.column_stack((cd_val_Cd0.flatten(),cd_val_Cd1.flatten(),cd_val_Cd5.flatten(),rmse_data.flatten(),meod_val.flatten()))
        print("Iteration ",str(cdind+1),"out of ", str(numcds.shape[0]))

# Bundle up final results from all fastmech runs. 
# Combine columns 1-5 from meas_wse csv file to measured and simulated WSEs data
meas_and_sim_wse2 = np.hstack((meas_wse[:,1:5],meas_and_sim_wse))
# save RMSE results file
np.savetxt(config.get('Params','rmse_file'), trmse, delimiter=',')
# Save measured vs simulated WSE data
np.savetxt(config.get('Params','meas_vs_sim_file'), meas_and_sim_wse2, delimiter=',')
print (rmse_data)
print (meas_and_sim_wse2)

# put results of the model rmse into a dataframe
df = pd.DataFrame(trmse,columns = ['Cd_Cd0','Cd_Cd1','Cd_Cd5','rmse','meon'])
# get the row of min rmse
best_mod = df.loc[df['rmse'].idxmin()]

# alternatively, get a certain run 
#best_mod = df.iloc[44,:]

# re-run the fastmech model with the best model parameters and save to working directory
Cd_Cd0 = best_mod['Cd_Cd0']
Cd_Cd1 = best_mod['Cd_Cd1']
Cd_Cd5 = best_mod['Cd_Cd5']
hdf5_file_name = "FM_Calib_" + str(Q) + "_CdCd0" + str(Cd_Cd0) + "_CdCd1" + str(Cd_Cd1) + "_CdCd5" + str(Cd_Cd5) + "_.CGNS"
copyfile(config.get('Params','base_file'),
         hdf5_file_name)
fastmech_BCs(hdf5_file_name, Q, H_DS, iniType, OneDCD, Cddistribution)
fastmech_params(hdf5_file_name, 1500)
fastmech_change_var_cd2(hdf5_file_name, newCd_Cd0 = Cd_Cd0,newCd_Cd1 = Cd_Cd1,newCd_Cd5 = Cd_Cd5)
for path in execute(["Fastmech.exe", hdf5_file_name]):
    meod = path
    try:
        meod = re.search('Discharge:(.+?)\n', meod).group(1)
    except AttributeError:
        meod = '-9999' # apply your error handling
    meod =  "".join(meod.split())
    meod = float(meod)
    
# Construct grid with solution values and create vtk modules or tools for extracting data by point
SGrid = vtk.vtkStructuredGrid()
create_vtk_structured_grid_post(SGrid, hdf5_file_name, xoffset, yoffset)
cellLocator2D = vtk.vtkCellLocator()
cellLocator2D.SetDataSet(SGrid)
cellLocator2D.BuildLocator()
WSE_2D = SGrid.GetPointData().GetScalars('WSE')
Depth_2D = SGrid.GetPointData().GetScalars('Depth')
Vel_2D = SGrid.GetPointData().GetScalars('Velocity')

# For the re-run model, eExtract WSEs for the given Q for the pixel center coordinate of the area of interest
# Create container for coordinate of the area of interest  
studyareaWSE_box = np.zeros(study_area_tbl.shape[0])
studyareaDepth_box = np.zeros(study_area_tbl.shape[0])
studyareaVel_box = np.zeros(study_area_tbl.shape[0])
for counter,line in  enumerate (study_area_tbl):
    point2D = [line[1]-xoffset, line[2]-yoffset, 0.0]
    pt1 = [line[1]-xoffset, line[2]-yoffset, 10.0]
    pt2 = [line[1]-xoffset, line[2]-yoffset, -10]
    idlist1 = vtk.vtkIdList()
    cellLocator2D.FindCellsAlongLine(pt1,pt2,0.0, idlist1)
    # retrieve simulated WSEs for study area
    try:
        WSEout = getCellValue(SGrid, point2D, idlist1.GetId(0), WSE_2D)
    except ValueError:
        WSEout = -9999
    studyareaWSE_box[counter] = WSEout
    try:
        Depthout = getCellValue(SGrid, point2D, idlist1.GetId(0), Depth_2D)
    except ValueError:
        Depthout = -9999
    studyareaDepth_box[counter] = Depthout
    try:
        Velout = getCellValue(SGrid, point2D, idlist1.GetId(0), Vel_2D)
    except ValueError:
        Velout = -9999
    studyareaVel_box[counter] = Velout
# put model WSE results for entire grid as a new column in the study_area_tbl
study_area_tbl_out  = pd.DataFrame(study_area_tbl,columns = ['ID','X','Y','Z'])
study_area_tbl_out['WaterSurfaceElevation'] = studyareaWSE_box
study_area_tbl_out['Velocity'] = studyareaVel_box
study_area_tbl_out['Depth'] = studyareaDepth_box
# save table as csv, if desired
#    study_area_tbl_out.to_csv(re.sub("\\.","_",str(Q))  + '_' + config.get('Params','study_area_tbl_out'),index = False)

# calculate average depth of inundation for all wetted cells (Depth != 0)
study_area_tbl_out
Inundated_cells = study_area_tbl_out[study_area_tbl_out['Depth'] > 0]
Depth_col = Inundated_cells[['Depth']]
mean_depth = Depth_col.mean(axis = 0)[0]
Vel_col = Inundated_cells[['Velocity']]
mean_Vel = Vel_col.mean(axis = 0)[0]
print("Average Depth of inundation: " + str(mean_depth))
print("Average Velocity: " + str(mean_Vel))
# LEV = (0.01 - 0.001)x average depth x average velocity
print("LEV range = " + str(0.01*mean_depth*mean_Vel) +' to ' + str(0.001*mean_depth*mean_Vel))
end = datetime.datetime.now()

elapsed = end - start
print('Time to run:' + str(elapsed))



