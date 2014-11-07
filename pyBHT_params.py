import numpy as np

# input and output files for BHT data and model results
input_file = 'input/BHTinput.csv'
outputfilename = 'results/BHTout.csv'

# file containing thermal and porosity parameters for user-specified litholgoies
litho_file = 'input/litho_params.csv'

# file containing thermal parameters of drilling mud
mud_file = 'input/mud_params.csv'

# file containing thermal parameters of pore water
water_file = 'input/water_params.csv'

# create a figure of the model results for each BHT series:
make_figure = True

# False for forward model only, True for calibration of model params:
calibrate = True

# if True: calibrate drilling mud temperature
# False: calibrate only the formation temperature, 
#  use fixed mud temperature from input file 
calibrate_mud_temp = True

# minimum value of mud temperature
# use this to constrain the search algorithm
# and prevent unrealistically low values
minimum_mud_temp = 10.0

# select inverse model algorithm
# choices: ['leastsq','simplex']
calibration_method = 'simplex'

# mesh size:
nr = 500
# model grid cell size (m)
dr = 0.005

# domain parameter, p=1 for radial, p=2 for spherical, p=0 for 1D
p = 1

# initial model timestep (sec)
# deprecated, model calculates optimal timestep
dt = 5.0

# drilling mud circulation time (hrs)
# this value is used if no value is 
# specified in the input .csv file for a BHT series
circtime_est = 5.0

# set this to true to equalize temperatures in the borehole during each
# temperature measurement
borehole_mixing = False