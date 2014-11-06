"""
calculate formation temperature from time series of three or more 
bottom hole temperatures taken at a single depth

simulates the cooling as a result of drilling simulated using a fixed
(low) temperature in the borehole and the subsequent thermal recovery

forward solution calculated using a 2D numerical finite difference 
model of heat flow

inverse modeling to estimate formation and borehole temperature 
using downhill simplex algorithm provided by scipy. 

2011-2014
Elco Luijendijk <elco.luijendijk@gmail.com>
""" 

import math
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
import scipy.optimize

import pyBHT_params as params
import lib.pyBHTlib as pyBHTlib


print 'Calculation of formation temperatures from time series of ' \
      'bottom hole temperatures'
print 'using a 2-media (borehole and formation) numerical ' \
      'finite-difference heat flow model'
print '2011-2014'
print 'Elco Luijendijk <elco.luijendijk@gmail.com>'

# read input datafile
df_raw = pd.read_csv(params.input_file)
#df_raw.set_index('id')
n_data = len(df_raw)

dfl = pd.read_csv(params.litho_file)
dfm = pd.read_csv(params.mud_file)
dfw = pd.read_csv(params.water_file)

if params.calibrate_mud_temp is True:
    n_parameters = 2
else:
    n_parameters = 1

# add  litho params to dataframe
df = df_raw.merge(dfl)

# add mud params to dataframe
df = df.merge(dfm)

# add water params to dataframe
df = df.merge(dfw)

# calculate porosity
df['porosity_final'] = df['porosity']
ind = df['porosity'].isnull()
df['porosity_final'][ind] = (df['phi_0'][ind]
                             * np.exp(-df['compressibility'][ind]
                                      * df['depth'][ind]))

# calculate matrix and formation thermal conductivity
df['K_formation'] = ((df['K_water']**df['porosity_final'])
                     * (df['K_matrix']**(1.0 - df['porosity_final'])))
df['density_formation'] = \
    (df['porosity_final'] * df['density_water']
     + ((1.0 - df['porosity_final']) * df['density_matrix']))
df['spec_heat_formation'] = \
    (df['porosity_final'] * df['specific_heat_water']
     + ((1.0 - df['porosity_final']) * df['specific_heat_matrix']))

# use estimated circulation time from parameter file
# if not specified in BHT input file
df['circ_time'][df['circ_time'].isnull()] = params.circtime_est

df = df.set_index(['id'])

# save input data + porosity and thermal parameters
print 'saving input data, calculated porosity and thermal parameters to %s' \
    % params.outputfilename
df.to_csv(params.outputfilename)

for well_no, bht_id in enumerate(df.index):

    print '\n--------- run %s / %s,  well %s ------\n' \
          % (well_no + 1,
             n_data,
             df.ix[bht_id, 'well'])

    # create arrays of BHT and recovery times
    BHTs = np.zeros(df.ix[bht_id, 'N_BHTs'])
    recovery_times = np.zeros(BHTs.shape)
    for bht_number in xrange(1, df.ix[bht_id, 'N_BHTs'] + 1):
        BHTs[bht_number-1] = df.ix[bht_id, 'T_%i' % bht_number]
        recovery_times[bht_number-1] = df.ix[bht_id,
                                             'rec_time_%i' % bht_number]

    # construct calibration parameters:
    parameters = np.zeros(n_parameters, dtype=float)
    parameters[0] = df.ix[bht_id, 'T_formation_init']
    if params.calibrate_mud_temp is True:
        parameters[1] = df.ix[bht_id, 'T_mud_init']

    if params.calibrate is True:

        if params.calibration_method == 'leastsq':
            return_data = True
        else:
            return_data = False

        args = (params.nr,
                params.dr,
                params.dt,
                df.ix[bht_id, 'circ_time'] * 60.0 * 60.0,
                df.ix[bht_id, 'diameter'] / 2.0,
                df.ix[bht_id, 'K_formation'],
                df.ix[bht_id, 'K_mud'],
                df.ix[bht_id, 'spec_heat_formation'],
                df.ix[bht_id, 'specific_heat_mud'],
                df.ix[bht_id, 'density_formation'],
                df.ix[bht_id, 'density_mud'],
                BHTs,
                recovery_times * 60.0 * 60.0,
                df.ix[bht_id, 'T_formation_init'],
                df.ix[bht_id, 'T_mud_init'],
                return_data,
                params.minimum_mud_temp,
                params.p,
                params.borehole_mixing)

        if params.calibration_method == 'leastsq':

            ## perform least squares optimization,  test: 8 runs,  min=1.66
            results = scipy.optimize.leastsq(pyBHTlib.objective_function,
                                             parameters,
                                             args=args,
                                             ftol=0.01)
            parameters = results[0]

            if params.calibrate_mud_temp is True:
                df.ix[bht_id, 'calibrated_formation_temp'] = parameters[0]
                df.ix[bht_id, 'calibrated_mud_temp'] = parameters[1]
            else:
                df.ix[bht_id, 'calibrated_formation_temp'] = parameters

            print 'calibrated parameters: ', parameters

        elif params.calibration_method == 'simplex':
            # optimize using simplex algorithm,  test: 1.66,  33 runs
            results = scipy.optimize.fmin(pyBHTlib.objective_function,
                                          parameters,
                                          args=args,
                                          xtol=0.1, ftol=0.1)
            #parameters = results
            df.ix[bht_id, 'calibrated_formation_temp'] = results[0]
            if params.calibrate_mud_temp is True:
                df.ix[bht_id, 'calibrated_mud_temp'] = results[1]

            print 'calibrated parameters: ', results

    # set final parameters:
    parameters = np.zeros(n_parameters, dtype=float)
    parameters[0] = df.ix[bht_id, 'calibrated_formation_temp']
    if params.calibrate_mud_temp is True:
        parameters[1] = df.ix[bht_id, 'calibrated_mud_temp']

    print '-' * 7
    print 'starting final model run with parameters: ', parameters

    # run simulation with final optimized parameters:
    results = pyBHTlib.simulate_BHTs(parameters,
                                     df.ix[bht_id, 'T_mud_init'],
                                     params.nr,
                                     params.dr,
                                     params.dt,
                                     df.ix[bht_id, 'circ_time'] * 60.0 * 60.0,
                                     df.ix[bht_id, 'diameter'] / 2.0,
                                     df.ix[bht_id, 'K_formation'],
                                     df.ix[bht_id, 'K_mud'],
                                     df.ix[bht_id, 'spec_heat_formation'],
                                     df.ix[bht_id, 'specific_heat_mud'],
                                     df.ix[bht_id, 'density_formation'],
                                     df.ix[bht_id, 'density_mud'],
                                     BHTs,
                                     recovery_times * 60.0 * 60.0,
                                     params.p,
                                     params.borehole_mixing,
                                     make_figure=params.make_figure)

    if params.make_figure is True:
        BHTout, RMSE,  BHT_times,  BHT_curve,  r, Tplot = results
    else:
        BHTout, RMSE = results
    
    # calculate R^2
    ssxy = np.sum((BHTs - BHTs.mean()) * (BHTout - BHTout.mean()))
    ssxx = np.sum((BHTs - BHTs.mean())**2)
    ssyy = np.sum((BHTout - BHTout.mean())**2)
    
    df.ix[bht_id, 'R2'] = ssxy**2 / (ssxx * ssyy)
    df.ix[bht_id, 'RMSE'] = RMSE

    for i, BHT in enumerate(BHTout):
        df.ix[bht_id, 'T_calibrated_%i' % (i+1)] = BHT

    print '--\nOptimized formation and mud temperatures of well %s:'\
          % df.ix[bht_id, 'well']
    print df.ix[bht_id, ['calibrated_formation_temp',
                         'calibrated_formation_temp']]

    # create figure:
    if params.make_figure is True:

        width = 110.0 / 25.4
        height = width

        fig = pl.figure(figsize=(width, height))

        fig_params = {'axes.labelsize': 'x-small',
                      'legend.fontsize': 'xx-small',
                      'xtick.labelsize': 'x-small',
                      'ytick.labelsize': 'x-small'}

        pl.rcParams.update(fig_params)

        axs = [fig.add_subplot(2, 1, 1),
               fig.add_subplot(2, 1, 2)]

        minval = int(math.floor(Tplot.min() / 2.) * 2)
        maxval = int(math.ceil(Tplot.max() / 2.) * 2)
        contourInt = np.arange(minval, maxval, 2)
        
        degree_symbol = unichr(176)
        axistext = 'Temperature (%sC)' % degree_symbol

        # temperature fields:
        for bhtNo in xrange(df.ix[bht_id, 'N_BHTs']):

            label = 'T, %0.1f hrs' % recovery_times[bhtNo]
            axs[0].plot(r, Tplot[bhtNo], lw=1.0, label=label)

        # temperature curve
        axs[1].plot(BHT_times / (60.0*60.0), BHT_curve, color='black', lw=1.0)
        axs[1].scatter(recovery_times, BHTs, facecolor='gray',
                       edgecolor='black')
        axs[1].set_ylim(BHT_curve.min(), BHT_curve.max() * 1.2)

        # make fill area for borehole
        radius = df.ix[bht_id, 'diameter'] / 2.0
        ylim = axs[0].get_ylim()
        axs[0].fill((0, radius, radius, 0),
                    (ylim[0], ylim[0], ylim[1], ylim[1]),
                    color='grey',
                    zorder=0,
                    label='borehole')

        tekst = 'R2 = %0.2f\nRMSE = %0.1f %sC' \
                % (df.ix[bht_id, 'R2'],
                   df.ix[bht_id, 'RMSE'],
                   degree_symbol)
        axs[1].text(0.9, 0.1, tekst,
                    ha='right', va='bottom',
                    fontsize='xx-small',
                    transform=axs[1].transAxes)

        axs[0].set_xlabel('Distance (m)')
        axs[1].set_xlabel('Recovery time (hr)')
        axs[0].set_ylabel(axistext)
        axs[1].set_ylabel(axistext)

        axs[0].legend(loc='lower right', fontsize='xx-small')

        for ax in axs:
            ax.set_xlim(0, ax.get_xlim()[-1])
            ax.grid()

        title = 'Temperature during recovery\n%s, well %s, depth=%0.0f m' \
                % (bht_id, df.ix[bht_id, 'well'], df.ix[bht_id, 'depth'])

        axs[0].set_title(title, fontsize='x-small')

        fig.tight_layout()

        fn = '%s_%s_%0.0fm_BHT_model.png' % (bht_id,
                                             df.ix[bht_id, 'well'],
                                             df.ix[bht_id, 'depth'])
        figfile = os.path.join('fig', fn)
        pl.savefig(figfile, dpi=300)

        print '-'*10
        print 'Figure saved as %s' % figfile
        print '-'*10

    # save input data + porosity and thermal parameters
    print 'saving input data and model results to %s' \
        % params.outputfilename
    df.to_csv(params.outputfilename)

print 'done'