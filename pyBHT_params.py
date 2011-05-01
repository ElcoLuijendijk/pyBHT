import numpy as np

# input and output files:
fileName = 'input/BHTinput.csv'
outputfilename = 'results/BHTout.csv'

# create a figure of the model results for each BHT series:
makeFigure = False

# False for forward model only, True for calibration of model params:
calibrate = False

# if True: calibrate drilling mud temperature
# False: calibrate only the formation temperature, 
#  use fixed mud temperature from input file 
calibrateMudTemp = True

# parameter that determines temperature averaging in the borehole
# 0: no averaging
# 1: temperatures averaged before each measurement
stir = 1

# select inverse model algorithm
# choices: ['leastsq','simplex']
optMethod = 'simplex'

# mesh size:
nx, ny = 200, 200
# model grid cell size (m)
cellsize = 0.01
# initial model timestep (sec)
# deprecated, model calculates optimal timestep
timestep = 5.0 
# drilling mud circulation time (hrs)
# this value is used if no value is 
# specified in the input .csv file for a BHT series
circtime_est = 5.0

# thermal conductivity:
# clay, sand, sand/clay mix, limestone, organic, evaporite
conductivity = np.array([ 2.5 , 3.5 ,  3.0 , 3.2 ,  1.0 ,  5.5 ])
# thermal conductivity drilling mud:
KMud = 0.88
# water:
Kwater = 0.6

# heat capacity
cRock = 0.9e+3 
cWater = 0.42e+3
cMud = 0.88e+3

# density
rhoRock = 2650.0   
rhoWater = 1025.0   
rhoMud = 1100.0 

# porosity-depth eq. data
# phi values:
# clay, sand, sand/clay mix, limestone, organic, evaporite
Phi0 = np.asarray([ 0.63 , 0.49 , 0.55 , 0.70 , 0.05 , 0.01 ])
Phic = np.asarray([ 0.51e-3 , 0.27e-3 , 0.40e-3 , 0.71e-3 , 0 , 0 ])

