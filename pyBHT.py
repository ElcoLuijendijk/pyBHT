# pyBHT
# calculate formation temperatures from bottom hole temperature data 
# Copyright 2011, Elco Luijendijk


"""
calculate formation temperature from time series of three or more 
bottom hole temperatures taken at a single depth

simulates the cooling as a result of drilling simulated using a fixed
(low) temperature in the borehole and the subsequent thermal recovery

forward solution calculated using a 2D numerical finite difference 
model of heat flow

inverse modeling to estimate formation and borehole temperature 
using downhill simplex algorithm provided by scipy. 

Apr 2011 
Elco Luijendijk <elco.luijendijk@gmail.com>
""" 

import math,  sys,  csv, pdb

import numpy as np
import pylab as pl
from scipy import optimize

from pyBHT_params import *

sys.path.append("./lib/")
import pyBHTlib

print 'Calculation of formation temperatures from time series of bottom hole temperatures'
print 'using a 2-media (borehole and formation) numerical finite-difference model'
print 'april 2011'
print 'Elco Luijendijk <elco.luijendijk@gmail.com>'

                                                    
header, wellindex, inputArray =\
        pyBHTlib.readCSVArray_id(fileName, delimiter=",")

Nwells,  dataLength = np.shape(inputArray)

if calibrateMudTemp == True:
    Nparameters = 2
else:
    Nparameters = 1


#create array to store simuation results
outputArray = np.zeros((Nwells, dataLength+5), dtype=float)
outputArray[:, 0:6] = inputArray[:, 0:6]   # copy input array into output array

parameters = np.zeros((Nparameters), dtype=float)

for wellNo in xrange(Nwells):

    print '\n--------- run %s /%s,  well %s ------\n'\
            %(wellNo+1, Nwells, wellindex[wellNo])
    
    ###### set parameters
    inch = 0.02504
    radius = inputArray[wellNo, 2]*inch/2.0
    Nbhts = int(inputArray[wellNo, 4])
    # circulation time in hrs
    circtimeHr = float(inputArray[wellNo, 3]) 
    # use the estimate entered by the user if no data available
    if circtimeHr == 999 or circtimeHr == 0:  
        circtimeHr = circtime_est
    # no. of timesteps circulation time
    circtime = int(circtimeHr*60.0*60.0) 
    
    ### set thermal parameters
    # calculate rock thermal properties using porosity-depth equation
    lithology = inputArray[wellNo, 1]
    # lithology from 4th column in input file 
    # 0=clay, 1=sand, 2=carbonate, 3=salt, 4=coal
    depth = inputArray[wellNo,  0] # depth below surface in m
    
    Phi0_bht = Phi0[lithology]
    Phic_bht = Phic[lithology]
    phi = Phi0_bht * np.exp(-Phic_bht*depth)
    
    # thermal conductivity W m-1 K-1 or J m-1 K-1 s-1
    Kmatrix = conductivity[lithology] 
    KRock = (Kwater**phi)*(Kmatrix**(1-phi))
    rhoRock_bulk = rhoRock*phi+rhoWater*(1-phi)
    cRock_bulk = cRock*phi+cWater*(1-phi)
        
    # set initial temperature estimates   
    initialTemp = float(inputArray[wellNo,  5]) 
    mudTemp = float(inputArray[wellNo,  6])

    # set BHT array
    BHTArray = np.zeros((2,  Nbhts), dtype=float)
    for bhtNo in range(Nbhts):
        BHTArray[0, bhtNo] = inputArray[wellNo,  7+bhtNo]
        BHTArray[1, bhtNo] = inputArray[wellNo,  14+bhtNo]
    recoveryTimes = BHTArray[0] 
    
    parameters = np.zeros((Nparameters), dtype=float)
    formationTemp = initialTemp
    parameters[0] = initialTemp
    if calibrateMudTemp == True:
        parameters[1] = mudTemp

    if calibrate == True:
        if optMethod == 'leastsq':
                
                ## perform least squares optimization,  test: 8 runs,  min=1.66
                returnData = True
                results = optimize.leastsq(pyBHTlib.residualFunc, parameters,
                                        args=(nx, ny, cellsize, timestep,
                                        circtime, radius, KRock, KMud,
                                        cRock_bulk, cMud, rhoRock_bulk, 
                                        rhoMud, stir, BHTArray,
                                        returnData), ftol=0.01)  
                parameters = results[0]
                if calibrateMudTemp == True:
                    formationTemp = parameters[0]
                    mudTemp = parameters[1]    
                else:
                    formationTemp = parameters
                
        elif optMethod == 'simplex':
            # optimize using simplex algorithm,  test: 1.66,  33 runs
            returnData = False
            results = optimize.fmin(pyBHTlib.residualFunc, parameters, 
                                    args=(nx, ny, cellsize, timestep,
                                    circtime, radius, KRock, KMud,
                                    cRock_bulk, cMud, rhoRock_bulk, rhoMud,
                                    stir, BHTArray, returnData),
                                    xtol=0.1, ftol=0.1)
            #parameters = results
            formationTemp = results[0]
            if calibrateMudTemp == True:
                mudTemp = results[1]    

    # set final parameters:
    parameters = np.zeros((Nparameters), dtype=float)
    parameters[0] = formationTemp
    if calibrateMudTemp == True:
        parameters[1] = mudTemp
    
    # run simulation with final optimized parameters:
    results = pyBHTlib.BHTcalcFunc(parameters, mudTemp, nx, ny,
                                    cellsize, timestep, circtime,
                                    radius, 
                                    KRock, KMud, cRock, cMud,
                                    rhoRock, rhoMud,
                                    stir, BHTArray,
                                    makeFigure=makeFigure)
    if makeFigure == True:
        BHTout, RMSE,  BHTtimes,  BHTcurve,  Tplot = results  
    else:
        BHTout, RMSE = results
    

    print '--\nOptimized formation and mud temperatures of well %s:'\
            %wellindex[wellNo]
    print ' %0.1f,  %0.1f\n--' %(formationTemp, mudTemp)

    # store results in an array:
    outputArray[wellNo,  11:(11+Nbhts)] = BHTout[2, :]
    outputArray[wellNo,  6] = KRock / (cRock_bulk*rhoRock_bulk)
    outputArray[wellNo,  7] = KMud / (cMud*rhoMud)
    outputArray[wellNo,  8] = formationTemp
    outputArray[wellNo,  9] = mudTemp
    outputArray[wellNo,  10] = RMSE

    # create figure:
    if makeFigure == True:
        
        Nrows = int(math.ceil((Nbhts+1)/2.0))
        
        pyBHTlib.initFigure(vert_size = Nrows*60.0)
        pl.subplots_adjust(wspace=0.25,  hspace=0.4)
    
        minval = int(math.floor(Tplot.min()/2.)*2)
        maxval = int(math.ceil(Tplot.max()/2.)*2)
        contourInt = np.arange(minval, maxval, 2)
        
        degree_symbol = unichr(176)
        axistext='Temperature (%sC)'%(degree_symbol)
        
        # temperature fields:
        for bhtNo in xrange(Nbhts):
            pl.subplot(Nrows, 2, bhtNo+1)
            
            im = pl.imshow(Tplot[:, :, bhtNo+1], vmin=minval,
                            vmax=maxval, cmap=pl.get_cmap('hot'))
            cn = pl.contour(Tplot[:, :, bhtNo+1], contourInt,
                            colors='gray')
            pyBHTlib.plotBoreholeRadius(radius, cellsize)
            pl.xlim(0, nx) ; pl.ylim(0, ny)
            
            titletxt = ['(A)', '(B)', '(C)', '(D)', '(E)', '(F)',
                        '(G)'][bhtNo]
            titletxt += ' %0.1f hrs recovery time'\
                %(recoveryTimes[bhtNo]/60.0)
            pl.title(titletxt)
            pl.xlabel('Distance (m)') ; pl.ylabel('Distance (m)')

        # temperature curve
        pl.subplot(Nrows, 2, bhtNo+2)
        pl.plot(BHTtimes/(60.0*60.0), BHTcurve, color='black', lw=1.0)
        pl.scatter(BHTArray[0]/60.0, BHTArray[1], facecolor='gray',
                    edgecolor='black')
        pl.xlabel('Time (hr)') ; pl.ylabel(axistext)
    
        titletxt = ['(A)', '(B)', '(C)', '(D)', '(E)', '(F)',
                    '(G)', '(H)'][bhtNo+1]
        titletxt += ' Borehole temperature'
        pl.title(titletxt)
        
        pl.subplots_adjust(hspace=0.5, wspace=0.3, bottom=0.17)
        
        axcb = pl.axes((0.1,0.07,0.4,0.02),frameon=False)
        pl.xticks([]) ; pl.yticks([])
        cb = pl.colorbar(im, cax=axcb,
                            orientation='horizontal')
        cb.add_lines(cn)
        cb.set_label(axistext, fontsize='xx-small')
            
    
        figfile = 'fig/'+wellindex[wellNo]+'_BHTmodel.png'
        pl.savefig(figfile,  dpi=300)
        print '-'*10
        print 'Figure saved as %s' %figfile
        print '-'*10
        
# save results:
header_output = header[0:7]+['diffusivity_rock', 'diffusivity_mud', 
    'calibrated formation temperature', 'calibrated mud temperature', 
    'RMSE', 'simulated BHT 1', '-2', '-3', '-4', '-5', '-6']
pyBHTlib.saveCSVArray_id(outputfilename, header_output, wellindex,
                            outputArray)

print '%i BHT recovery calculations done' %(wellNo+1)
print 'results saved as csv file:\n%s' %(outputfilename)
