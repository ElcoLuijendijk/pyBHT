# pyBHT
# calculate formation temperatures from bottom hole temperature data 
# Copyright 2011, Elco Luijendijk

import math, sys, csv, pdb

import numpy as np
import pylab as pl
from scipy import optimize

sys.path.append("./lib/")
import heatflow, heatflow_v2

print heatflow_v2.__doc__

def initFigure(hor_size=190.0, vert_size=-1, textsize='x-small'):
    
    golden_ratio = (1.0+math.sqrt(5))/2.0
        
    if vert_size == -1:
        vert_size = hor_size/golden_ratio
    elif vert_size == 1:
        vert_size = hor_size*golden_ratio
        
    # initialize figure
    pl.figure(figsize=(hor_size/25.4, vert_size/25.4))  # set figure size to a4 
    
    # set default parameters for figure
    if textsize == 'xx-small':
        textsize_s = 'xx-small'
        textsize_l = 'x-small'
        textsize_leg = 'xx-small'
    elif textsize == 'x-small':
        textsize_s = 'x-small'
        textsize_l = 'small'
        textsize_leg = 'xx-small'
    elif textsize == 'small':
        textsize_s = 'small'
        textsize_l = 'medium'
        textsize_leg = 'xx-small'
    params = {'axes.labelsize': textsize_s, 'text.fontsize': textsize_l\
    , 'legend.fontsize': textsize_leg, 'axes.titlesize' : textsize_l\
    , 'xtick.labelsize': textsize_s, 'ytick.labelsize': textsize_s}
    
    pl.rcParams.update(params)
    
    return


def readCSVArray_id(fileName, delimiter = '\t'):
    
    # function to read a csv file
    fin = open(fileName, 'r')
    #fileContents  = fin.read()
    f = csv.reader(fin, delimiter=delimiter)
    well_id = [] ; data = [] ; rowNo = 0 
    for row in f:
        # Save header row.
        if rowNo == 0:
            header = row
        else:
            data.append(row[1:])        
            well_id.append(row[0])
        rowNo += 1
    fin.close()
    # convert to array
    dataArray = np.asarray(data).astype(float)
    return header, well_id, dataArray

def saveCSVArray_id(fileName, header, wellIndex, data, delimiter=','):
    
    fout = open(fileName, 'w')
    f = csv.writer(fout, delimiter=delimiter)
    f.writerow(header)
    Nrows, Ncols = np.shape(data)
    for row in xrange(Nrows):
        f.writerow([wellIndex[row]]+data[row].tolist())
    
    fout.close()
    
    return


def BHTcalcFunc(parameters, mudTemp, nx, ny, cellsize, timestep,
        circtime, radius, KRock, KMud, cRock, cMud, rhoRock, rhoMud, 
        stir, BHTArray, makeFigure=False, useFipy=False, debug=False):

    """
    function to calculate BHTs
    uses FD solver function for temperature recovery model
    """
    
    # initialize thermal parameter grids
    T = np.zeros((nx, ny), dtype=float)
    KGrid = np.zeros((nx, ny), dtype=float)
    cGrid = np.zeros((nx, ny), dtype=float)
    rhoGrid = np.zeros((nx, ny), dtype=float)
    bound = np.zeros((nx, ny), dtype=float)
    Nrows, Nbhts = np.shape(BHTArray)

    if makeFigure ==True:
        BHTcurve = np.array([])
        BHTtimes = np.array([])
        Tplot = np.zeros((nx, ny, Nbhts+1))
    
    # set up thermal parameter grids
    print 'start setting up temperature and thermal parameters of\
            %sx%s grid' %(nx, ny)
    for x in range(0, nx):
        for y in range(0, ny):
            # calculate distance from 0, to see whether in borehole or formation
            distance = (math.sqrt((x*cellsize)**2+(y*cellsize)**2))
            if distance > radius:
                KGrid[x, y] = KRock
                cGrid[x, y] = cRock
                rhoGrid[x, y] = rhoRock
                T[x, y] = parameters[0]  
            else:
                KGrid[x, y] = KMud
                cGrid[x, y] = cMud
                rhoGrid[x, y] = rhoMud
                try:
                    T[x, y] = parameters[1]
                except:
                    T[x, y] = mudTemp
                bound[x, y] = 1    # boundary condition
    borehole = bound.copy()
    
    # calculate max timestep size
    diffusivity = KGrid / (cGrid*rhoGrid)
    maxTimeStep = cellsize**2 / (2*diffusivity.max())
    timeStep_adj = maxTimeStep / 2.0
    Nsteps = int(circtime/timeStep_adj)
    print '\ttime steps: %i, time step size = %0.2e sec' %(Nsteps, timeStep_adj)
    
    ########################################
    # simulate temperatures during drilling:
    ########################################
    print '\tsimulate T, drilling mud circulation'
    if makeFigure == True:
        BHTline = np.zeros(Nsteps)
        T, BHTline = heatflow_v2.heatflow_v2(T, BHTline, 
                            KGrid, cGrid, rhoGrid, cellsize,
                            timeStep_adj, bound, borehole,
                            Nsteps, nx, ny)
        # store data for figure:
        #BHTcurve = np.concatenate((BHTcurve, BHTline))
        #try:
        #    startTime = BHTtimes[-1] 
        #except:
        #    startTime = 0
        #BHTtimes = (startTime+np.linspace(timeStep_adj,
        #                    timeStep_adj*Nsteps, Nsteps))
            
        Tplot[:, :, 0] = T.copy()
    else:
        T = heatflow.heatflow(T, KGrid, cGrid, rhoGrid, cellsize, Nsteps,
                        timeStep_adj,
                        bound, nx, ny)
        
        if debug == True:
            V = np.arange(int(T.min()), T.max()+2.0, 2.0)
            pl.clf()
            pl.subplot(3,2,1)
            cf = pl.contourf(T,V)
            pl.colorbar()
     
            
    ####################################################
    # simulate temperature recovery after drilling
    ####################################################

    # remove temperature boundary condition:
    bound[:, :] = 0
    # create array to store simulated BHTs
    BHTout = np.zeros((Nrows+1, Nbhts), dtype=float)    
    # copy BHT input array into output Array
    BHTout[:-1, :] = BHTArray      
    totalTime = 0 ; sqerror = 0
    print '\tsimulate T, recovery'
    for i in xrange(Nbhts):
        recovery = BHTArray[0, i] * 60.0
        timeleft = recovery-totalTime
        totalTime = recovery
        # set no of timesteps:
        Nsteps = int(timeleft/timeStep_adj)
        # simulate temperature evolution:
        
        if makeFigure ==True:
            BHTline = np.zeros(Nsteps)
            T, BHTline = heatflow_v2.heatflow_v2(T, BHTline, 
                                KGrid, cGrid, rhoGrid, cellsize,
                                timeStep_adj, bound, borehole,
                                Nsteps, nx, ny)
            # store temperature data for figure:
            try:
                startTime = BHTtimes[-1] 
            except:
                startTime = 0
            BHTcurve = np.concatenate((BHTcurve, BHTline))
            BHTtimes = np.concatenate((BHTtimes,
                        (startTime+np.linspace(timeStep_adj,
                            timeStep_adj*Nsteps, Nsteps))))
            Tplot[:, :, 1+i] = T.copy()
        
        else:
            T = heatflow.heatflow(T, KGrid, cGrid, rhoGrid, cellsize,
                                Nsteps, timeStep_adj, bound, nx, ny)
            
        if debug == True:
            pl.subplot(3, 2, 3+i)
            pl.contourf(T, V)
            pl.colorbar()
        
        #calculate avg., max T in borehole section
        BHTavg = ((borehole*T).sum())/(borehole.sum())
            
        # simulate convection in borehole:
        # all temperatures in borehole are averaged
        # before recording temperature
        if stir == 1 and bound.max() == 0:
            T = -((borehole-1)*T) + BHTavg * borehole                   
                
        print '\tBHT %i of %i, simulated T:  %0.2f, obs. BHT: %0.2f'\
                %(i+1, Nbhts, BHTavg, BHTArray[1, i])
        
        # store simulated BHT in output array BHTout
        BHTout[2, i] = BHTavg 
        sqerror += (BHTavg-BHTArray[1, i])**2
    
    # save figure of model results:
    if debug == True:
        pl.savefig('t_field.png')
        
    # calculate RMSE of observed and simulated BHTs:
    RMSE = math.sqrt(sqerror/Nbhts)
    
    if makeFigure == True:
        return BHTout, RMSE, BHTtimes, BHTcurve, Tplot
    else:
        return BHTout, RMSE


def testpoint(radius, cellsize, xt, yt):
    xt_=xt-(0.5*cellsize)
    yt_=yt-(0.5*cellsize)
    distance=math.sqrt(xt_**2+yt_**2)
    if distance<=radius:
        argument=True
    else:
        argument=False
    #print 'radius= %s, dist= %s, %s' %(radius, distance, argument)
    
    return argument


def plotBoreholeRadius(radius, cellsize):
    x = math.ceil(radius/cellsize)*cellsize ; y=0
    xpoints=[] ; ypoints = []
    count=0
    while x>=0:
        xpoints.append(x/cellsize) ; ypoints.append(y/cellsize)
        yt = y + cellsize 
        if testpoint(radius, cellsize, x, yt) == True:
            y = yt
        else:
            x = x - cellsize 
    xpoints.append(x/cellsize) ; ypoints.append(y/cellsize)
    pl.plot(np.asarray(xpoints), np.asarray(ypoints), color='k', 
            linewidth=2)        
    
    return


def residualFunc(parameters, nx, ny, cellsize, timestep,
    circtime, radius, KRock, KMud, cRock, cMud, rhoRock, rhoMud, stir,
    BHTArray, returnData):
    
    if len(parameters) == 2:
        print 'calibration params: T formation: %0.2f, T mud: %0.2f'\
            %(parameters[0], parameters[1])
    elif len(parameters) == 1:
        print 'calibration params: T formation: %0.2f' %(parameters[0])
    
    initialTemp=parameters[0]
    mudTemp=parameters[1]
    BHTout, RMSE = BHTcalcFunc(parameters, mudTemp, nx, ny,
                                cellsize, timestep, circtime, radius,
                                KRock, KMud, cRock, cMud, rhoRock,
                                rhoMud, stir, BHTArray)

    #residuals=zeros(shape(BHTout), dtype=float)
    residuals = BHTout[2, :]-BHTArray[1, :]
        
    print 'RMSE of observed and simulated T: %0.2f \n ------' %RMSE
    
    if returnData == True:
        return residuals
    else:
        return RMSE
    
