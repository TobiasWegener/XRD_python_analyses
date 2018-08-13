# -*- coding: utf-8 -*-
'''
Created on Wed Aug 08 2018
@author: Tobias Wegener

Modul with the different functions needed for the XRD analyses jupyter notebook
it contains:

Functions:
__________

json_data_read_reduce: 
function to read in XRD data as jason file with can be downloaded from materialsproject

peak: 
function that calculates the lattice constance for bcc metals and alloys from the peak position, here for W and Cr

relevant: 
filters the peak position which are condsideren in the callulation of the peak fuction above

irrelevant: 
filters all peaks, which were not considered above for the oxids peaks

ox_index: 
searchsorted based functions that compares two arrays, in a uncertanty interval and gives back the index of values from the array one which fit to array two with in the uncertanty. This makes it possible to filter only the fitting peaks to the litereature data.

plane_count: 
gives information about the number of alloy peaks and their distribution. This helps to find the optimal number of bins for the sorting of alloy peaks in to different concentrations.

'''

version = '1.5'

import numpy as np



def json_data_read_reduce(file_path, start = 20, stop = 120, min_amp=20):
    '''read in json files from "https://www.materialsproject.org" for xrd data,
    and transform them to one numpy array of 2theta data for comparison.
    
    Parameters
    ----------
    file_path:  file path string with ending e.g.: "Y2WO6_mp-565796_xrd_Cu.json"
    start:      start of relevant positions in degree, default 20
    stop:       stop of relevant position, default 120
    min_amp:    minmal amplitude of values considered
    
    Returns
    -------
    two_theta value array
    
    Example
    ------
    Y6WO12 = json_data_read_reduce('daten_bank/Y6WO12_mp-19005_xrd_Cu.json')
    print 'Y6WO12', len(Y6WO12), Y6WO12
    >>> Y6WO12 6 [28.42636114 29.20471251 33.62374459 47.78369588 48.78884306 57.43435056]
    '''
    #import dependeces
    import json
    import numpy as np
    
    #open file and assine to d
    with open(file_path) as json_data:
        d = json.load(json_data)
        
    #produce array with the same length as data
    amplitude = np.zeros(len(d[u'pattern'][:]))
    two_theta = np.zeros(len(d[u'pattern'][:]))
    #loop through the length of data and assine (fill) the arrays with the data
    for i in range(len(d[u'pattern'][:])):
        amplitude[i] = d[u'pattern'][i][0]
        two_theta[i] = d[u'pattern'][i][2]
    
    #filter only revant, data is of interest
    start_stop = (start <= two_theta) & (two_theta <= stop)
    #print two_theta[start_stop]
    
    #filter only data with relevant amplitude here, default higher or equal 20%
    minimal_amp = (amplitude > min_amp)
    #print amplitude[minimal_amp]
        
    #return the filtered array of the json file which fullfills minimal amplitude and range in 2Theta
    return two_theta[start_stop & minimal_amp]



def peak(x,lambd=1.5418):
    '''Calculates the lattice constance for bcc metals and alloys
    and only takes one plane after the other of the XRD 
    to fit the constance. To avoid mix up of planes. 
    The used order and the  given boundarys are (110) with 40.23 <x) & (x<44.41), 
    (200)(58. <x) & (x< 65.6), (211)(73.17 <x) & (x< 81.78),(220) (86.9 <x) & (x< 98.04) 
    for W-Cr alloys. 
    
    Parameters
    ----------
    x = fitted position in 2Theta (°)
    lambd = wavelength of the Kathode metal in angstrom default is Cu. If other  source is used
    or if Ka2 is stripe, give the correct value
    
    Retun
    ------
    Tuple yf np.arrays of fitted constant'''
    
    a, b, c, d, e, f, g, h = 40.23, 44.41, 58.1, 64.6, 73.17, 81.78, 86.9, 98.04
    x = np.array(x)
    x = x[(a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)]
    #Lattice constant a for the (110) plane
    a_110 = np.array(lambd/np.sin(x[(a <x) & (x < b)]*np.pi/180/2)/2**0.5)
    a_200 = np.array(lambd/np.sin(x[(c <x) & (x < d)]*np.pi/180/2))
    #appen calculated a's for panes (110) and (200)
    a_12 = np.append(a_110, a_200)
    a_211 = np.array(0.5*lambd/np.sin(x[(e < x) & (x < f)]*np.pi/180/2)*6**0.5)
    a_220 = np.array(0.5*lambd/np.sin(x[(g <x) & (x< h)]*np.pi/180/2)*8**0.5)
    a_34 = np.append(a_211, a_220)
    #appen all in one np.array
    a = np.append(a_12,a_34)
    return a



def relevant(x):
    '''Calculates the peak position which are really considered in the calculation above return list of values, Further more it gives back the list of the index
    
    Returns
    _______
    r the array of filtered values which are relevant
    and y an array of binary values to use the filter on other arrays with the same length'''
    
    a, b, c, d, e, f, g, h = 40.23, 44.41, 58.1, 64.6, 73.17, 81.78, 86.9, 98.04
    x = np.array(x)
    r =  x[(a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)]
    y =   [(a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)]
    return r, y



def irrelevant(x):
    '''Calculates the peak position which should not be considered in the calculate above
    
    Returns
    -------
    List of values, Further more it gives back the list of the index'''
    a, b, c, d, e, f, g, h = 40.23, 44.41, 58.1, 64.6, 73.17, 81.78, 86.9, 98.04
    x = np.array(x)
    #r =  x[(a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)]
    y =   (a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)
    y = np.logical_not(y)
    return y



#searchsorted numpy based fitting ---> very fast and precise /prone if values overlap in the uncertanty range
def ox_index(a, b, uncertanty = 0.07, sort = False, unique = False, reduce_to_relevant = False):
    r'''Searchsorted based functions, that compares two arrays, 
    in an uncertanty interval (default [-0.07, 0.07]) and gives back the index of values 
    from the array one which fit to array two. 
    This makes it possible to filter, only the fitting peaks to the litereature data. 
    Base on the numpy.searchsorted function
    (https://docs.scipy.org/doc/numpy/reference/generated/numpy.searchsorted.html).
    
    Parameters
    ----------
    a : array of values (here peak postisions)
    b : array of values to compare to (here, litreture, angels to expect)
    uncertanty : expected uncertanty, where the data should fit in [-lower,upper]
    
    Returns
    -------
    index: Index of vlaues in the array a vlues which fit b array
    
    Notes
    -----
    gives back measurd values, which fit to the database (in the given uncertanty range) both array have to be sorted!!
    
    Examples
    --------
    measured_values = np.array([12.01,13.001,14.001,15.001,15.6,15.9999])
    database = np.array([12.,13.,14.,15.,16.])
    print measured_values[Ox_index(measured_values,database)]
    >>> [ 12.01    13.001   14.001   15.001   15.9999]
    
    #It also works in both directions:
    a = np.array([1.01,2.01,2.5,3.01,3.999,4.0001,6,7,9,89])
    b = np.array([1.,2.,3.,4.,6])
    print 'a',a[Ox_index(a,b)]
    print 'b',b[Ox_index(b,a)]
    print np.abs(a[Ox_index(a,b)]-b[Ox_index(b,a)]).mean()
    >>>a [ 1.01    2.01    3.01    3.999   4.0001  6.    ]
    >>>b [ 1.  2.  3.  4.  4.  6.]
    >>>0.00518333333333
    '''
    #import dependences
    import numpy as np
    
    #check which data type and correct to np.array if needed
    if type(a) != "<type 'numpy.ndarray'>":
        a = np.array(a)
    if type(b) != "<type 'numpy.ndarray'>":
        b = np.array(b)
    
    # if wanted sort input array a other wise search sort evetnualy fails
    if sort == True:
        a = np.sort(a)
        b = np.sort(b)
    
    # optinal, reduce input arrays to unique values:
    if unique == True:
        a = np.unique(a)
        b = np.unique(a)
    
    # checks if distance in array b is smaller then uncertanty to avoid overlapping errors
    if uncertanty > np.ediff1d(b).min():
        print 'Uncertanty {}, smallest distance in the array is {:.3f}-> updated uncertanty!'.format(uncertanty ,np.ediff1d(b).min())
        uncertanty = np.ediff1d(b).min()
    
    # optional, reduces the input arrays, so that the longer array a values, 
    #fall in the range of the short array b and visa versa
    if reduce_to_relevant == True:
        a = a[((b.min()-uncertanty)<a) & (a<(b.max()+uncertanty))]
        b = b[(a.min()<(b-uncertanty)) & ((b+uncertanty)<a.max())]
    
    def intersection(a,b):
        '''main function: finding indexes which are the intersection of the 
        lower and higher uncertanty range of a in b'''
        x = np.searchsorted(a,b-uncertanty,side='left')
        y = np.searchsorted(a,b+uncertanty,side='right')-1
        index = x[x==y]
        return np.unique(index) 
    
    index = intersection(a,b)
    #retuns the index of the fitting parts of the array with further information:
    if len(index) == 0:
        #print ('No match found for uncertanty {:.3f}'.format(uncertanty))
        return index
    elif len(intersection(b,a)) == len(index):
        if np.abs(b[intersection(b,a)]-a[index]).mean()<(uncertanty-0.014):#(len(intersection(a,b))*uncertanty/10)):
            print('{} peaks mean deviation from data is {:.3f}'.format(len(index),np.abs(b[intersection(b,a)]-a[index]).mean()))
            return index
        print ('{} peaks do unknown uncertantiy, may high strain, check?'.format(len(intersection(a,b))))
        return index
    else:
        print ('{} peaks unknown deviation but in the uncertantiy {:.3f}'.format(len(intersection(a,b)),uncertanty))
        return index



def plane_count(x):
    '''
    Follwing should be true:
    1. The number of peaks for one concentration can not be higher than the number of Fitted planes 
        => active_planes >= fequanzy musst be lower or equal.-> So bin number has to be increased in this case
    2. The number of concentrations can't be lower, than the highest number of peaks fitted in one plane, 
       depends on the lower boundary of peak fitts required
    This helps to find a reasonable amount of bins and the highest frequancy the physic allows.
    
    Returns
    -------
    active_planes :   number of fitted planes as 
    minimal__phases : minimal number of phases to be seperated, based of the highest peak count in one plane 
    fitted_planes :   count per plane    
    
    Example:
    --------
    x = np.array([42,45,50,59,73.5,87,88,100])
    print 'highest possible number of frequancy (number of planes):', plane_count(x)[0] #condition 1.
    print 'minimal number of phases:', plane_count(x)[1] #condition 2.
    >>> highest possible number of frequancy (number of planes): 4 
    >>> minimal number of phases: 2
    '''
    #number of fitted planes
    x = np.array(x)
    a, b, c, d, e, f, g, h = 40.23, 44.41, 58.1, 64.6, 73.17, 81.78, 86.9, 98.04
    fitted_planes = np.array([len(x[(a <x) & (x< b)]),len(x[(c <x) & (x< d)]),len(x[(e <x) & (x< f)]),len(x[(g <x) & (x< h)])])
    #condition 1. number of active planes
    active_planes = (fitted_planes >= 1).sum()    
    #condition 2. number of minimal number of concentration
    minimal_phases = fitted_planes.max() 
        
    return active_planes, minimal_phases, fitted_planes
