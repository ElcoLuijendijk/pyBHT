*****
PyBHT
*****

Description
===========
PyBHT is a model code that calculates formation temperature from time
series of three or more bottom hole temperatures recorded at the same
depth. The model simulates cooling as a result of drilling and the
subsequent thermal recovery using an explicit finite difference 
solution of the radial heat flow equation.
Formation and borehole temperature are calibrated using the downhill 
simplex algorithm provided by Scipy.


Usage
=====
- Use a text editor to adjust the location of the BHT input file and 
  the model parameters in the model parameter file "PyBHT_params.py"
- Construct a .csv file containing bottom hole temperature input data.
  An example file is located in: /input/BHT.csv
- The thermal parameters of the formation rocks surrounding the borehole,
  pore water, and drilling mud are located in separate .csv files in the folder
  /input.
- Run the model code::
        
        python PyBHT.py

- A .csv file that contains the model results, including the calibrated
  formation temperatures is stored in the folder /results/

- A figure of the model results for each BHT series is created automatically
  and stored in the folder /fig/


Dependencies
============
PyBHT requires the following python packages:

* Python 2.x: http://www.python.org/

* NumPy: http://www.scipy.org/NumPy

* Matplotlib: http://matplotlib.sourceforge.net/

* SciPy: http://www.scipy.org/


Reference
=========

Please cite the following article if you publish work that uses PyBHT:

Luijendijk, E., M. Ter Voorde, R.T. Van Balen, H. Verweij, E. Simmelink. (2011)
Thermal state of the Roer Valley Graben, part of the European Cenozoic Rift System.
Basin Research, 23(1), 65-82.
DOI: 10.1111/j.1365-2117.2010.00466.x

A bibtex file of this citation can be found in the PyBHT folder

You can find a copy of the paper at: 
http://dx.doi.org/10.1111/j.1365-2117.2010.00466.x. 
Contact me by email if you do not have access to this journal.


License
=======

PyBHT is distributed under the GNU General Public License, version 3:
http://www.gnu.org/copyleft/gpl.html
A copy of this license is distributed along with PyBHT, see: gpl-3.0.txt

Elco Luijendijk <elco.luijendijk at gmail.com>
Nov 2014


