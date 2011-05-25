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

import math,  sys,  csv, pdb, os

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

# compile fortran heat flow modules:
if len(sys.argv)>1 and sys.argv[1] == 'compile':
    os.chdir('./lib')
    os.system('f2py -m -c heatflow heatflow.f')
    os.system('f2py -m -c heatflow_v2 heatflow_v2.f')
    os.chdir('./..')

                                                    
header, wellindex, inputArray =\
        pyBHTlib.readCSVArray_id(fileName, delimiter=",")

Nwells,  dataLength = np.shape(inputArray)

if calibrateMudTemp == True:
    Nparameters = 2
else:
    Nparameters = 1


#create array to store simuation results
Ninp = 8
outputArray = np.zeros((Nwells, Ninp + 13 ), dtype=float)
outputArray[:, :Ninp] = inputArray[:, :Ninp]   # copy input array into output array

parameters = np.zeros((Nparameters), dtype=float)

for wellNo in xrange(Nwells):

    print '\n--------- run %s /%s,  well %s ------\n'\
            %(wellNo+1, Nwells, wellindex[wellNo])
    
    ###### set parameters
    radius = inputArray[wellNo, 5]
    Nbhts = int(inputArray[wellNo, 7])
    # circulation time in hrs
    circtimeHr = float(inputArray[wellNo, 6]) 
    # use the estimate entered by the user if no data available
    if circtimeHr >= 99999 or circtimeHr == 0:  
        circtimeHr = circtime_est
    # circulation time in seconds
    circtime = int(circtimeHr*60.0*60.0) 
    
    ### set thermal parameters
    # calculate rock thermal properties using porosity-depth equation
    lithology = inputArray[wellNo, 1]
    # lithology from 4th column in input file 
    # 0=clay, 1=sand, 2=carbonate, 3=salt, 4=coal
    depth = inputArray[wellNo,  0] # depth below surface in m

    # porosity
    phi = inputArray[wellNo, 2]    
    if phi == 0 or phi >= 99999:
        Phi0_bht = Phi0[lithology]
        Phic_bht = Phic[lithology]
        phi = Phi0_bht * np.exp(-Phic_bht*depth)
    
    diffusivity_rock = inputArray[wellNo, 3]
    
    if diffusivity_rock == 0 or diffusivity_rock >= 99999:
        # thermal conductivity W m-1 K-1 or J m-1 K-1 s-1
        Kmatrix = conductivity[lithology] 
        KRock = (Kwater**phi)*(Kmatrix**(1-phi))
        rhoRock_bulk = rhoRock*phi+rhoWater*(1-phi)
        cRock_bulk = cRock*phi+cWater*(1-phi)

    else:
        # if diffusivity is specified:
        # calculate thermal conductivity from diffusivity
        # TODO: change input params of FD solver so that only
        # diffusivity param is required as input
        rhoRock_bulk = rhoRock*phi+rhoWater*(1-phi)
        cRock_bulk = cRock*phi+cWater*(1-phi)
        KRock = diffusivity_rock / (cRock_bulk * rhoRock_bulk)
    
    diffusivity_mud = inputArray[wellNo, 4] 
    if diffusivity_mud != 0 and diffusivity_mud < 99999:
        print '\tusing specified diffusivity of drilling mud'
        # thermal conductivity of drilling mud from diffusivity
        KMud2 = diffusivity_mud / (cMud * rhoMud )
    else:
        # use standard thermal conductivity defined in parameter file
        KMud2 = KMud
    
    # set initial temperature estimates   
    initialTemp = float(inputArray[wellNo,  8]) 
    mudTemp = float(inputArray[wellNo,  9])

    # set BHT array
    BHTs = np.zeros((Nbhts), dtype=float)
    recoveryTimes = np.zeros((Nbhts), dtype=float)
    for bhtNo in range(Nbhts):
        BHTs[bhtNo] = inputArray[wellNo,  17+bhtNo]
        recoveryTimes[bhtNo] = inputArray[wellNo,  10+bhtNo]
    
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
                                        circtime, radius, KRock, KMud2,
                                        cRock_bulk, cMud, rhoRock_bulk, 
                                        rhoMud, stir, BHTs,
                                        recoveryTimes,
                                        mudTemp, minimumMudTemp,
                                        returnData),
                                        ftol=0.01)  
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
                                    circtime, radius, KRock, KMud2,
                                    cRock_bulk, cMud, rhoRock_bulk,
                                    rhoMud, stir, BHTs, recoveryTimes,
                                    mudTemp, minimumMudTemp,
                                    returnData),
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
                                    KRock, KMud2, cRock, cMud,
                                    rhoRock, rhoMud,
                                    stir, BHTs, recoveryTimes,
                                    makeFigure=makeFigure)
    if makeFigure == True:
        BHTout, RMSE,  BHTtimes,  BHTcurve,  Tplot = results  
    else:
        BHTout, RMSE = results
    
    # calculate R^2
    ssxy = ((BHTs - BHTs.mean()) * (BHTout - BHTout.mean())).sum()
    ssxx = (((BHTs - BHTs.mean())**2).sum())
    ssyy = (((BHTout - BHTout.mean())**2).sum())
    
    R2 = ssxy**2 /(ssxx * ssyy)

    print '--\nOptimized formation and mud temperatures of well %s:'\
            %wellindex[wellNo]
    print ' %0.1f,  %0.1f\n--' %(formationTemp, mudTemp)

    # store results in an array:
    outputArray[wellNo,  (Ninp+6):(Ninp+6+Nbhts)] = BHTout[:]
    outputArray[wellNo,  Ninp] = KRock / (cRock_bulk*rhoRock_bulk)
    outputArray[wellNo,  Ninp+1] = KMud2 / (cMud*rhoMud)
    outputArray[wellNo,  Ninp+2] = formationTemp
    outputArray[wellNo,  Ninp+3] = mudTemp
    outputArray[wellNo,  Ninp+4] = RMSE
    outputArray[wellNo,  Ninp+5] = R2

    # create figure:
    if makeFigure == True:
        
        Nrows = int(math.ceil((Nbhts+1)/2.0))
        
        pyBHTlib.initFigure(vert_size = Nrows*60.0)
        pl.subplots_adjust(wspace=0.25,  hspace=0.4)
    
        minval = int(math.floor(Tplot[:,:,1:].min()/2.)*2)
        maxval = int(math.ceil(Tplot[:,:,1:].max()/2.)*2)
        contourInt = np.arange(minval, maxval, 2)
        
        degree_symbol = unichr(176)
        axistext='Temperature (%sC)'%(degree_symbol)
        
        # temperature fields:
        for bhtNo in xrange(Nbhts):
            pl.subplot(Nrows, 2, bhtNo+1)
            
            im = pl.imshow(Tplot[:, :, bhtNo+1], vmin=minval,
                            vmax=maxval, cmap=pl.get_cmap('hot'))
            cn = pl.contour(Tplot[:, :, bhtNo+1], 
                            colors='gray')
            pyBHTlib.plotBoreholeRadius(radius, cellsize)
            pl.xlim(0, nx) ; pl.ylim(0, ny)
            
            titletxt = ['(A)', '(B)', '(C)', '(D)', '(E)', '(F)',
                        '(G)'][bhtNo]
            titletxt += ' %0.1f hrs recovery time'\
                %(recoveryTimes[bhtNo])
            pl.title(titletxt)
            pl.xlabel('Distance (m)') ; pl.ylabel('Distance (m)')

        # temperature curve
        axc = pl.subplot(Nrows, 2, bhtNo+2)
        pl.plot(BHTtimes/(60.0*60.0), BHTcurve, color='black', lw=1.0)
        pl.scatter(recoveryTimes, BHTs, facecolor='gray',
                    edgecolor='black')
        pl.ylim(BHTcurve.min(),BHTcurve.max()*1.2)
        axc.yaxis.grid(True)
        axc.xaxis.grid(True)
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
header_output = header[:Ninp+1]+['diffusivity_rock', 'diffusivity_mud', 
    'T_calibrated', 'Tmud_calibrated', 
    'RMSE', 'R2', 'BHTsim_1', 'BHTsim_2', 'BHTsim_3', 'BHTsim_4', 'BHTsim_5',
    'BHTsim_6', 'BHTsim_7']
pyBHTlib.saveCSVArray_id(outputfilename, header_output, wellindex,
                            outputArray)

print '%i BHT recovery calculations done' %(wellNo+1)
print 'results saved as csv file:\n%s' %(outputfilename)
