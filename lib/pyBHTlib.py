# pyBHT
# calculate formation temperatures from bottom hole temperature data 
# Copyright 2011-2014, Elco Luijendijk

import math
import pdb
import numpy as np


def simulate_BHTs(parameters, mud_temp, nr, dr, dt,
                  circtime, radius,
                  K_formation, K_mud,
                  c_formation, c_mud,
                  rho_formation, rho_mud,
                  BHTs, recovery_times,
                  p,
                  borehole_mixing,
                  make_figure=False,
                  debug=False):

    """
    Simulate borehole cooling and thermal recovery
    """
    
    # initialize thermal parameter grids
    r = np.arange(0, dr*nr, dr) + dr
    r_mid = (r[1:] + r[:-1]) / 2.0
    T = np.zeros(nr)
    K_array = np.zeros(nr-1)
    c_array = np.zeros(nr)
    rho_array = np.zeros(nr)

    n_bhts = np.shape(BHTs)[0]

    # set up thermal parameters in FD mesh
    print '\tsetting up initial T and parameters'
    borehole = r <= radius
    formation = borehole == False

    K_array[borehole[1:]] = K_mud
    c_array[borehole] = c_mud
    rho_array[borehole] = rho_mud

    K_array[formation[1:]] = K_formation
    c_array[formation] = c_formation
    rho_array[formation] = rho_formation

    T[formation] = parameters[0]

    if len(parameters) > 1:
        mud_temp = parameters[1]

    T[borehole] = mud_temp

    # generate q array
    q = np.zeros((nr + 1))

    # calculate max timestep size
    c_mid = (c_array[1:] + c_array[:-1]) / 2.0
    rho_mid = (rho_array[1:] + rho_array[:-1]) / 2.0
    diffusivity = K_array / (c_mid * rho_mid)
    max_timestep = dr**2 / (2 * diffusivity.max())
    dt_adj = max_timestep / 2.0
    n_steps = int(circtime / dt_adj)
    print '\ttime steps: %i, time step size = %0.2e sec' \
          % (n_steps, dt_adj)
    
    ########################################
    # simulate temperatures during drilling:
    ########################################
    print '\tsimulate T, drilling mud circulation'

    for step in xrange(n_steps):
        T = radial_explicit_heat_flow(T, q, r, r_mid, dr,
                                      K_array, rho_array, c_array,
                                      dt_adj, p)

        T[borehole] = mud_temp

    # calculate avg., max T in borehole section
    BHTavg = T[borehole].mean()
        
    print '\t + %0.2f hr, T= %0.2f' \
          % (n_steps * dt_adj / (60.0*60.0), BHTavg)

    ####################################################
    # simulate temperature recovery after drilling
    ####################################################

    # create array to store simulated BHTs
    BHTout = np.zeros(n_bhts,
                      dtype=float)

    total_time = 0
    sqerror = 0

    if make_figure is True:
        T_plot = np.zeros((n_bhts, T.shape[0]))
        BHT_temp = []
        BHT_times = []

    print '\tsimulate T, recovery'
    for i in xrange(n_bhts):
        recovery = recovery_times[i]
        runtime = recovery - total_time
        total_time += runtime
        # set no of timesteps:
        n_steps = int(runtime / dt_adj)
        # simulate temperature evolution:

        if make_figure is True:

            BHT_temp_rec = np.zeros(n_steps)
            BHT_time_rec = np.arange(n_steps) * dt_adj + total_time - runtime

            for step in xrange(n_steps):
                T = radial_explicit_heat_flow(T, q, r, r_mid, dr,
                                              K_array, rho_array, c_array,
                                              dt_adj, p)

                BHT_temp_rec[step] = T[borehole].mean()
                T_plot[i] = T.copy()

            BHT_temp.append(BHT_temp_rec)
            BHT_times.append(BHT_time_rec)

        else:
            for step in xrange(n_steps):
                T = radial_explicit_heat_flow(T, q, r, r_mid, dr,
                                              K_array, rho_array, c_array,
                                              dt_adj, p)

        if borehole_mixing is True:
            T[borehole] = T[borehole].mean()

        #calculate avg., max T in borehole section
        BHTavg = T[borehole].mean()
            
        print '\t+%0.2f hr,BHT %i/%i, sim. T:%0.2f, obs. BHT:%0.2f' \
              % (n_steps*dt_adj / (60.0*60.0), i+1, n_bhts,
                 BHTavg, BHTs[i])
        
        # store simulated BHT in output array BHTout
        BHTout[i] = BHTavg 
        sqerror += (BHTavg - BHTs[i])**2

    # merge BHT time series
    if make_figure is True:
        BHT_temp_array = np.concatenate(BHT_temp)
        BHT_times_array = np.concatenate(BHT_times)

    # calculate RMSE of observed and simulated BHTs:
    RMSE = math.sqrt(sqerror / n_bhts)

    if debug is True:
        pdb.set_trace()

    if make_figure is True:
        return BHTout, RMSE, BHT_times_array, BHT_temp_array, r, T_plot
    else:
        return BHTout, RMSE


def objective_function(parameters, nr, cellsize, timestep,
                       circtime, radius,
                       K_formation, K_mud,
                       c_formation, c_mud,
                       rho_formation, rho_mud,
                       BHTs,
                       recovery_times,
                       formation_temperature, mud_temperature,
                       return_data,
                       minimum_mud_temp,
                       p,
                       borehole_mixing):

    """
    Simulate borehole cooling and recovery and return difference
     between observed and simulated bottom hole temperatures
    """

    if len(parameters) == 2:
        print 'calibration params: T formation: %0.2f, T mud: %0.2f'\
            % (parameters[0], parameters[1])
        mud_temperature = parameters[1]

        if minimum_mud_temp is not None and mud_temperature < minimum_mud_temp:
            mud_temperature = minimum_mud_temp
            parameters[1] = minimum_mud_temp

    elif len(parameters) == 1:
        print 'calibration params: T formation: %0.2f' %(parameters[0])
    
    initial_temp = parameters[0]
    
    BHTout, RMSE = simulate_BHTs(parameters, mud_temperature, nr,
                                 cellsize, timestep, circtime, radius,
                                 K_formation, K_mud,
                                 c_formation, c_mud,
                                 rho_formation, rho_mud,
                                 BHTs, recovery_times,
                                 p, borehole_mixing)

    #residuals=zeros(shape(BHTout), dtype=float)
    residuals = BHTout - BHTs
        
    print 'RMSE of observed and simulated T: %0.2f \n ------' % RMSE
    
    if return_data is True:
        return residuals
    else:
        return RMSE


def radial_explicit_heat_flow(T, q, r, r_mid, dr, K, rho, cp, dt, p):

    """
    explicit finite difference solution of radial heat flow eq.
    """

    q[1:-1] = -K * r_mid**p * np.diff(T) / dr
    dT = (q[:-1] - q[1:]) / dr * (dt / (r**p * rho * cp))
    T += dT

    return T


def explicit_heat_flow_1d(T, q, dx, K, rho, cp, dt):

    """
    explicit finite difference solution of 1D heat flow eq.
    """

    q[1:-1] = -K * np.diff(T) / dx
    dT = (q[:-1] - q[1:]) / dx * (dt / (rho * cp))
    T += dT

    return T


def explicit_heat_flow_2d(T, qh, qv, dx, dy, Kh, Kv, rho, cp, dt):

    """
    explicit finite difference solution of 2D heat flow eq.
    """

    qh[1:-1, :] = -Kh * (T[1:, :] - T[:-1, :]) / dx
    qv[:, 1:-1] = -Kv * (T[:, 1:] - T[:, :-1]) / dy

    dTh = (qh[:-1, :] - qh[1:, :]) / dx * (dt / (rho * cp))
    dTv = (qv[:, :-1] - qv[:, 1:]) / dy * (dt / (rho * cp))

    dT = dTh + dTv

    T += dT

    return T