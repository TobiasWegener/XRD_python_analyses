## XRD_python_analyses
*Tool to analys XRD data for BCC metals and alloys. 
It fits the the peaks with help of the Package 'Peakutiles' and than calculates the lattice constant of the alloy phases and the formed oxides.*
## Description
### Two basic functions:
1. From the fited peaks, it caclulates the values for Cr concentration by the Vegards approximation, filters them by bins, and than shows these concentrations in the Plot legend.
#### Example for BCC phase fit:
![Alt text](https://github.com/TobiasWegener/XRD_python_analyses/blob/master/06_XRD_newLowPO2_bulk.svg "BCC phase fit")

2. It compares the rest of the peaks which are not detected as Metal/Alloy to a dictionary of different possible oxides, marks them in the plot and gives them in the legend. Here, a new very fast numpy compare algorithm was developted for this purpose, which makes it possible to compare to arrays of different length on approximate values.
#### Example for oxide and BCC phase fit:
![Alt text](https://github.com/TobiasWegener/XRD_python_analyses/blob/master/06_XRD_newBulk_rough_smoth.svg "BCC phase fit")

## Requirements
This repository contains Jupyter notebook (`.ipynb`) files. All of the scripts were written using Python 2.7. 
The main libraries used are [numpy](http://www.numpy.org/), pyplot from [matplotlib](https://matplotlib.org/index.html), [scipy.signal](https://docs.scipy.org/doc/scipy/reference/signal.html) with function 'savgol_filter' for smoothing, [peakutils](https://pypi.python.org/pypi/PeakUtils) for the fit, [sys](https://docs.python.org/2/library/sys.html) ,and [os](https://docs.python.org/2/library/os.html) for shell coments.

## Literature data
For the metal and oxides the data of ([ICDD-Database (2014)](www.icdd.com) spefically: PDF2014 PDF-2) is used,
with the exception of Y$_6$WO$_12$ which is taken from: Persson, K. (2014). Materials data on y6wo12 (sg:148) by [materials
project](https://www.materialsproject.org/materials/mp-19005/).

## Installation and use
Only download the (`.ipynb`) file, put an file in to a folder for example '/XRD', than add a folder of samples with the name discribing the exmperiment e.g. 'Bulk' in this folder all plots are saved. In this folder add a subfolder with the name 'data', here you can put the ('.xy') data files to investigate. The full path than looks like ('/XRD/Bulk/data).

## Coment and Feedback

At the moment this is optimized for the BCC W-Cr system but an addaption for other system is possible.
