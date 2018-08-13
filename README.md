## XRD_python_analyses
*Tool to analyze XRD data for BCC metals and alloys.
It fits the peaks with help of the Package 'Peakutiles' and then calculates the lattice constant of the alloy phases and the formed oxides.*
## Description
### Two basic functions:
1. From the fitted peaks, the values for the Cr concentration are calculated by the Vegard's approximation, the concentrations are filtered by bins to seperate the different BCC phases, and then the concentrations (weight%) are displayed in the Plot legend under the BCC title.

#### Example for BCC phase fit:
![Alt text](https://github.com/TobiasWegener/XRD_python_analyses/blob/master/06_XRD_newLowPO2_bulk.svg "BCC phase fit")

2. The rest of the peaks, which were not detected as Metal/Alloy, are compared to a dictionary of different possible oxides, then these are marked in the plot and given in the legend. Herein, a new very fast numpy compare algorithm, based on [numpy.searchsorted](https://docs.scipy.org/doc/numpy/reference/generated/numpy.searchsorted.html), was developted for this purpose making it possible to compare to arrays of different length on approximate values. It is available in the XRD_Analyses_modul.py and is called ox_index.
#### Example for oxide and BCC phase fit:
![Alt text](https://github.com/TobiasWegener/XRD_python_analyses/blob/master/06_XRD_newBulk_rough_smoth.svg "BCC phase fit")

## Requirements
This repository contains a Jupyter notebook ('.ipynb') file and the XRD_Analyses_Modul ('.py') file. All the scripts were written using Python 2.7, portation to Python 3.X is possible.
The main libraries used are [numpy](http://www.numpy.org/), pyplot from [matplotlib](https://matplotlib.org/index.html), [scipy.signal](https://docs.scipy.org/doc/scipy/reference/signal.html) with function 'savgol_filter' for smoothing, [peakutils](https://pypi.python.org/pypi/PeakUtils) for the fit, [sys](https://docs.python.org/2/library/sys.html), and [os](https://docs.python.org/2/library/os.html) for shell comments.

## Literature data
For the metal and oxides the data of ([ICDD-Database (2014)](www.icdd.com) spefically: PDF2014 PDF-2) is used,
except Y$_6$WO$_{12}$ and W$_{8}$O$_{21}$, which are taken from: Persson, K. (2014). Materials data on y6wo12 (sg:148) by [materials
project](https://www.materialsproject.org/materials/mp-19005/).

## Installation and use
Only download the ('.ipynb') and the XRD_Analyses_Modul files, put both files in to a folder for example '/XRD', then add a folder of samples with the name describing the experiment e.g. 'Bulk' in this folder all plots are saved. In this folder ('Bulk') create a subfolder with the name '/data', here you can put the ('.xy') data files to investigate. The full path then looks like ('/XRD/Bulk/data).

## Comment and Feedback

At the moment this is optimized for the BCC W-Cr system but an adaption for other systems is possible if the [Vegard's law](https://en.wikipedia.org/wiki/Vegard%27s_law) applies for the system.
Distributed under the MIT license. See ``LICENSE.txt`` for more information.
