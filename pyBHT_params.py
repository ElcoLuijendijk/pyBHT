import numpy as np

# input and output files:
fileName = 'input/BHTinput.csv'
outputfilename = 'results/BHTout.csv'

# create a figure of the model results for each BHT series:
makeFigure = True

# False for forward model only, True for calibration of model params:
calibrate = True

# if True: calibrate drilling mud temperature
# False: calibrate only the formation temperature, 
#  use fixed mud temperature from input file 
calibrateMudTemp = True

# minimum value of mud temperature
# use this to constain the search algorithm 
# and prevent unrealistic values
minimumMudTemp = 10.0

# parameter that determines temperature averaging in the borehole
# 0: no averaging
# 1: temperatures averaged before each measurement
stir = 1

# select inverse model algorithm
# choices: ['leastsq','simplex']
optMethod = 'simplex'

# mesh size:
nx, ny = 100, 100
# model grid cell size (m)
cellsize = 0.02
# initial model timestep (sec)
# deprecated, model calculates optimal timestep
timestep = 5.0 
# drilling mud circulation time (hrs)
# this value is used if no value is 
# specified in the input .csv file for a BHT series
circtime_est = 5.0

# thermal conductivity:
# sand, sand/clay mix, clay, limestone, organic, evaporite
conductivity = np.array([ 3.5 , 3.0 ,  2.5 , 3.2 ,  1.0 ,  5.5 ])
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
# sand, sand/clay mix, clay, limestone, organic, evaporite
Phi0 = np.asarray([  0.49 , 0.55 , 0.63, 0.70 , 0.05 , 0.01 ])
Phic = np.asarray([  0.27e-3 , 0.40e-3 , 0.51e-3, 0.71e-3 , 0 , 0 ])


