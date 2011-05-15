*****
PyBHT
*****

Description
===========
PyBHT is a model code that calculates formation temperature from time
series of three or more bottom hole temperatures recorded at the same
depth. The model simulates cooling as a result of drilling and the
subsequent thermal recovery using an explicit finite difference 
solver of the heat flow equation in a 2D model space.
Formation and borehole temperature are calculated using the downhill 
simplex algorithm provided by Scipy.


Usage
=====
- Use a text editor to adjust the location of the BHT input file and 
  the model parameters in the model parameter file "PyBHT_params.py"
- Construct a file containing the BHT series input data.
  An example file is located in: /input/BHT.csv
- If you are running the model code for the first time:

    - Compile the fortran heat flow modules using f2py::
        
        cd lib
        f2py -m -c heatflow heatflow.f
        f2py -m -c heatflow_v2 heatflow_v2.f
        
    - Or instruct PyBHT to compile the fortran modules by passing a 'compile'
      argument to the python script::
        
        python PyBHT.py compile

- Run the model code:: 
        
        python PyBHT.py

- A .csv file that contains the model results, including the calibrated 
  formation temperatures is stored in the folder /results/

- A figure of the model results is created automatically and stored in the folder /fig/


Dependencies
============
PyBHT requires the following packages to be installed on the system:


* Python: http://www.python.org/

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

You can find a copy at: 
http://dx.doi.org/10.1111/j.1365-2117.2010.00466.x. 
Contact me by email if you do not have access to this journal.


License
=======

PyBHT is distributed under the European Union Public Licence V. 1.1
A copy of this license is distributed along with PyBHT, see
license_EN.pdf and license_NL.pdf



Elco Luijendijk <elco.luijendijk at gmail.com>
April 2011 


