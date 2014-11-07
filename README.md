# PyBHT

## Description

PyBHT is a model code that calculates formation temperature from time
series of three or more bottom hole temperatures (BHT) recorded at the same
depth. The model simulates cooling as a result of drilling and the
subsequent thermal recovery using an explicit finite difference 
solution of the radial heat flow equation.
Formation and borehole temperature are calibrated using the downhill 
simplex algorithm provided by Scipy.


## Installation

- Clone pyBHT or download and unpack the zip file at the right hand side
- Navigate to the pyBHT folder and run the model code from a terminal, IDLE
  or your favorite python editor::

        python PyBHT.py

- The model will now calibrate the formation and mud temperature for 2 example
  BHT datasets from the Roer Valley Graben in the Netherlands


## Usage

- Use a text editor or excel to construct a .csv file containing bottom hole
  temperature input data. An example file is located in: `/input/BHTinput.csv`
- Use a text editor to adjust the default model parameters or the names of the 
  BHT input files in the model parameter file `PyBHT_params.py`
- The thermal parameters of the formation rocks surrounding the borehole,
  pore water, and drilling mud are located in the files: 
  `input/litho_params.csv`, `input/water_params.csv` and
  `input/mud_params.csv`. Adjust these if needed.
- For each BHT series you can specify a lithology type, pore water type and
  drilling mud type in the file `/input/BHTinput.csv`. The model will look for 
  the lithology, water type and mud type in the thermal parameter .csv files 
  and assign thermal parameters accordingly.
- The model calculates the bulk thermal conductivity, density and specific heat
  of the formation using data on porosity and the thermal parameters of the
  pore water and formation specified in `input/litho_params.csv` and
  input/water_params.csv. Porosity can be assigned directly for each BHT
  series in the file `input/BHTinput.csv`. 
- If porosity is not specified in `input/BHTinput.csv` (ie. the porosity column
  is left empty), porosity for each BHT series will be calculated using an 
  exponential porosity-depth equation: 
  phi = phi0 * exp(-compressibility * depth)
- The phi_0 and compressibility parameters for each lithology are specified in 
  the file `input/litho_params.csv`
- After running the model a .csv file that contains the model results,
  including the calibrated formation temperatures is stored in the folder
  `results/`. See `results/BHTout.csv` for an example of the output
- A figure of the model results for each BHT series is created automatically
  and stored in the folder `fig/`


## Dependencies

PyBHT requires the following python packages:

* Python 2.x: http://www.python.org/

* NumPy: http://www.scipy.org/NumPy

* Matplotlib: http://matplotlib.sourceforge.net/

* SciPy: http://www.scipy.org/

* Pandas: http://pandas.pydata.org/


## Reference

Please cite the following article if you publish work that uses PyBHT:

Luijendijk, E., M. Ter Voorde, R.T. Van Balen, H. Verweij, E. Simmelink. (2011)
Thermal state of the Roer Valley Graben, part of the European Cenozoic Rift System.
Basin Research, 23(1), 65-82. DOI: 10.1111/j.1365-2117.2010.00466.x

A bibtex file of this citation can be found in the PyBHT folder

You can find a copy of the paper at: 
http://dx.doi.org/10.1111/j.1365-2117.2010.00466.x. 

Contact me by email if you do not have access to this journal.


## License

PyBHT is distributed under the GNU General Public License, version 3:
http://www.gnu.org/copyleft/gpl.html

A copy of this license is distributed along with PyBHT, see: gpl-3.0.txt

Elco Luijendijk <elco.luijendijk at gmail.com>

November 2014


