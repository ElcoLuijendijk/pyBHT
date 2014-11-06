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

# minimumm temperature for drilling mud
# lower bound during claibration, use this to avoid unrealistically low mud T
minimum_mud_temp = 10.0

# set this to true to equalize temperatures in the borehole during each
# temperature measurement
borehole_mixing = False

# thermal conductivity for lithology numbers specified in BHT input file
# default lithologies:
# sand, sand/clay mix, clay, limestone, organic material, evaporite
#K = np.array([3.5, 3.0, 2.5, 3.2, 1.0, 5.5])
# thermal conductivity drilling mud:
#K_mud = np.array([0.88])
# water:
#K_water = 0.6

# heat capacity
#spec_heat = np.array([900.0, 900.0, 900.0, 900.0, 900.0, 900.0])
#spec_heat_water = 4181.3
#spec_heat_mud = np.array([880.0])

# density
#density = np.array([2650.0, 2650.0, 2650.0, 2650.0, 2650.0, 2650.0])
#density_water = 1025.0
#density_mud = np.array([1100.0])

# porosity-depth eq. data
# used in combination with athy's law to estimate porosity if not specified in
# BHT input file
# Athy's law: porosity = phi_0 * exp(-compressibility * depth)
# default lihtologies:
# [sand, sand/clay mix, clay, limestone, organic, evaporite]
#phi_0 = np.array([0.49, 0.55, 0.63, 0.70, 0.05, 0.01])
#compressibility = np.array([0.27e-3, 0.40e-3, 0.51e-3, 0.71e-3, 0, 0])


