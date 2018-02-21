# XRD_python_analyses
Tool to analys XRD data for BCC metals and alloys. 
It fits the the peaks with help of the Package 'Peakutiles' and than calculates the lattice constant of the alloy phases and the formed oxides.
## Two basic functions:
1. From the fited peaks, it caclulates the values for Cr concentration by the Vegards approximation, filters them by bins, and than shows these concentrations in the Plot legend.
2. It compares the rest of the peaks which are not detected as Metal/Alloy to a dictionary of different possible oxides, marks them in the plot and gives them in the legend. Here, a new very fast numpy compare algorithm was developted for this purpose, which makes it possible to compare to arrays of different length on approximate values.

At the moment this is optimized for W-Cr system but an addaption for other system ist possible.

