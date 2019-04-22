# Code for running a series of fastmech models, each with a different discharge and downstream boundary condition.
# Code reads in a csv file (QHDS.csv) that contains columns for Q and the downstream WSE. The code copies the base file (Case1.cgn)
# that is PRECONFIGURED by the user to grid conditions mapped (ie roughness) and calculation conditions set (roughness distribution variable by node)
# The code would need to be modified slightly to run with roughness distribution constant.
# Tool is developed by Julian Scott, but is strongly dependent on code written by Rich McDonald, especially the VTK, h5py utilization, and fastmech execution.


# 04/10/19



from __future__ import print_function  # Only Python 2.x
import numpy as np
from shutil import copyfile
import sys
import h5py
import vtk
import subprocess
import os
from itertools import count
import configparser
import re
import glob
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

#hdf_file = hdf5_file_name

def fastmech_BCs(hdf_file, iniType, OneDCD, Cddistribution):
    file = h5py.File(hdf_file, 'r+')
    group2 = file['/iRIC/CalculationConditions/FM_HydAttCDType/Value']
    dset2 = group2[u' data']
    dset2[0] = Cddistribution
    group4 = file['/iRIC/CalculationConditions/FM_HydAttWSType/Value']
    dset4 = group4[u' data']
    dset4[0] = iniType
    group7 = file['/iRIC/CalculationConditions/FM_HydAttWS1DCD/Value']
    dset7 = group7[u' data']
    dset7[0] = OneDCD
    file.close()
    


def fastmech_change_cd(hdf_file, newCd):
    # hdf5_file_name = r'F:\Kootenai Project\USACE\Braided\Case11_tmp.cgn'
    # r+ adds read/write permisions to file
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttCD/Value']
    dset = group[u' data']
    # print dset[0]
    dset[0] = newCd
    # print dset[0]
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
            dset2[index] = newWSE_0
    file.close()
    
def fastmech_change_H_DS(hdf_file, newH_DS):
    # hdf5_file_name = r'F:\Kootenai Project\USACE\Braided\Case11_tmp.cgn'
    # r+ adds read/write permisions to file
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttWS/Value']
    dset = group[u' data']
    dset[0] = newH_DS
    file.close()

def fastmech_change_Q(hdf_file, newQ):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_HydAttQ/Value']
    dset = group[u' data']
    dset[0] = newQ
    file.close()

def fastmech_params(hdf_file, Itermax, endLEV):
    file = h5py.File(hdf_file, 'r+')
    group = file['/iRIC/CalculationConditions/FM_EndLEV/Value']
    dset = group[u' data']
    dset[0] = endLEV
    group = file['/iRIC/CalculationConditions/FM_SolAttItm/Value']
    dset = group[u' data']
    dset[0] = Itermax
    file.close()


def add_fastmech_solver_to_path(solverPath):
    # print(os.environ['PATH'])
    os.environ['PATH'] += solverPath
    # print("\n")
    # print('new path')
    # print(os.environ['PATH'])

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
    ibc_grp = file['iRIC/iRICZone/FlowSolution1/IBC']
    velx_grp = file['iRIC/iRICZone/FlowSolution1/VelocityX']
    vely_grp = file['iRIC/iRICZone/FlowSolution1/VelocityY']

    xcoord_data = xcoord_grp[u' data']
    ycoord_data = ycoord_grp[u' data']
    wse_data = wse_grp[u' data']
    topo_data = topo_grp[u' data']
    ibc_data = ibc_grp[u' data']
    velx_data = velx_grp[u' data']
    vely_data = vely_grp[u' data']
    # SGrid = vtk.vtkStructuredGrid()
    ny, nx, = xcoord_data.shape
    # print(ny, nx)
    sgrid.SetDimensions(nx, ny, 1)
    points = vtk.vtkPoints()
    wseVal = vtk.vtkFloatArray()
    wseVal.SetNumberOfComponents(1)
    ibcVal = vtk.vtkIntArray()
    ibcVal.SetNumberOfComponents(1)
    velVal = vtk.vtkFloatArray()
    velVal.SetNumberOfComponents(1)
    for j in range(ny):
        for i in range(nx):
            points.InsertNextPoint(xcoord_data[j, i] - xoffset, ycoord_data[j, i] - yoffset, 0.0)
            wseVal.InsertNextValue(wse_data[j, i])
            ibcVal.InsertNextValue(ibc_data[j, i])
            velVal.InsertNextValue(np.sqrt(np.power(velx_data[j, i],2) + np.power(vely_data[j,i],2)))
        sgrid.SetPoints(points)

        sgrid.GetPointData().AddArray(wseVal)
        sgrid.GetPointData().AddArray(ibcVal)
        sgrid.GetPointData().AddArray(velVal)
    wseVal.SetName("WSE")
    ibcVal.SetName("IBC")
    velVal.SetName("Velocity")

# Set the config file 
# flow folder
q_folder = "q3_8_to_8_0"
setfile = "E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\demo_Directory\\Multiple_Discharge_Code\\" + q_folder + "\\Q_HDS_config.ini"
# Set up configuration file parser
config = configparser.ConfigParser()
# Read in the configuration file
config.read(setfile)


# GET XY points to extract WSE from each run
ObsWse_list = glob.glob("E:\\_DoD\\_Camp_Pendleton_Survey\\IRIC\_Modeling_dir\\_LowFlows_Model_v2\\Python_Directory\\demo_Directory\\Multiple_Discharge_Code\\" + q_folder + "\\ReadObsWSE_XYZ*.csv")

# =============================================================================
# Read in table of discharges and associated H_DS 
QHDS_table = np.genfromtxt(config.get('Params','QHDS_table'), delimiter=',', skip_header=0)
# Get number of Qs to run in batch
Q_count = sum(1 for row in QHDS_table) 
# =============================================================================

# Roughness distribution type; val = 1 for roughness variable by node 
Cddistribution = 1 

# Temporarily add fastmech solver to environment paths
add_fastmech_solver_to_path(config.get('Params','lib_path'))
add_fastmech_solver_to_path(config.get('Params','solver_path'))

# Set working directory from parameter in config file
os.chdir(config.get('Params','working_dir'))
# This is from Rich's original code; Leave set to 0.
xoffset = 0
yoffset = 0
# Initial condition type 
iniType = config.getfloat('Params','iniType')
# Drag coef for 1D Step backwater initial condition model
OneDCD = config.getfloat('Params','OneDCD')

# Solution iterations 
sol_iterations = 50

# Number of iterations to consider in rmse
numIters = 1
# Create output box (columns: Q, H_DS, Cd, RMSE, MEOD; rows = fastmech runs)
outputbox = np.zeros(shape=(Q_count,QHDS_table.shape[0]+2))
# get start time for timing batch process
start = datetime.datetime.now()

# For testing or stepping through for loop 
Qi = 0
row = 0
for Qi,row in enumerate(np.arange(0, Q_count, 1)):
    Q = QHDS_table[row,0]
    H_DS = QHDS_table[row,1]
    print("row: ",row)
    print("Q: ",Q)
    print("H_DS: ",H_DS)
    Q_match = "_" + re.sub("\\.","_",str(Q)) + ".csv"
    ObsQ_match = [s for s in ObsWse_list if Q_match in s]
    meas_wse = np.genfromtxt(ObsQ_match, delimiter=',', skip_header=1)
    nummeas = meas_wse.shape[0]
    meas_and_sim_wse = np.zeros(shape=(nummeas, numIters+1))
    hdf5_file_name = "FM_Calib_Q_" + re.sub("\\.","_",str(Q)) + ".CGNS"
    copyfile(config.get('Params','base_file'),hdf5_file_name)
    fastmech_change_H_DS(hdf5_file_name, H_DS)
    fastmech_change_Q(hdf5_file_name,Q)
    fastmech_BCs(hdf5_file_name,iniType, OneDCD, Cddistribution)
    fastmech_params(hdf5_file_name,sol_iterations, 0.0075)
    for path in execute(["Fastmech.exe", hdf5_file_name]):
#        print(path, end="")
        # Retrieve iteration output as a string
        meod = path
        # Extract from iteration output string the meon error on discharge;
        # If the solution is not resolved, return error value
        try:
            meod = re.search('Discharge:(.+?)\n', meod).group(1)
        except AttributeError:
            meod = '-9999' # apply your error handling
        # The final string in the solution iteration output contains the mean error on discharge.
        meod =  "".join(meod.split())
        meod = float(meod)
        
    # Construct grid with solution values and create vtk modules or tools for extracting data by point
    SGrid = vtk.vtkStructuredGrid()
    create_vtk_structured_grid(SGrid, hdf5_file_name, xoffset, yoffset)
    cellLocator2D = vtk.vtkCellLocator()
    cellLocator2D.SetDataSet(SGrid)
    cellLocator2D.BuildLocator()
    WSE_2D = SGrid.GetPointData().GetScalars('WSE')
    IBC_2D = SGrid.GetPointData().GetScalars('IBC')
    Velocity_2D = SGrid.GetPointData().GetScalars('Velocity')
    
    # Create container for simulated WSEs at the coordinates of the Obs WSEs    
    simwse = np.zeros(meas_wse.shape[0])
    # Create container for Obs WSEs
    measwse = np.zeros(meas_wse.shape[0])
    
    # For each ObWSE point, extract the simulated WSE 
    for counter, line in enumerate(meas_wse):
#        print(counter,line)
        point2D = [line[0]-xoffset, line[1]-yoffset, 0.0]
        pt1 = [line[0]-xoffset, line[1]-yoffset, 10.0]
        pt2 = [line[0]-xoffset, line[1]-yoffset, -10]
        idlist1 = vtk.vtkIdList()
        cellLocator2D.FindCellsAlongLine(pt1,pt2,0.0, idlist1)
        cellid = idlist1.GetId(0)
        # Get WSE cell value
        tmpwse = getCellValue(SGrid, point2D, cellid, WSE_2D)
        # populate the simulated WSE container at position 'counter'
        simwse[counter] = tmpwse
        # populate the Observed WSE container at position 'counter'
        measwse[counter] = line[2]
    rmse_i = rmse(simwse, measwse)
    outputbox[Qi,0] = Q
    outputbox[Qi,1] = H_DS
    outputbox[Qi,2] = rmse_i
    outputbox[Qi,3] = meod
    print("Discharge:",Q)
    print("Mean Error On Discharge:",meod)
    print("RMSE: ",rmse_i)
   
# put results into a dataframe and save to working directory
df = pd.DataFrame(outputbox,columns = ['Q','H_DS','RMSE','MEOD'])
df.to_csv(config.get('Params','outputbox'),index = False)

end = datetime.datetime.now()
print("Start time: " +str(start) + "   End time: " + str(end))

