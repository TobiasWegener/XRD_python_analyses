
# coding: utf-8

# ## Import Libraries

# In[61]:


get_ipython().magic(u'matplotlib inline')
import numpy as np
from scipy.signal import savgol_filter
import peakutils # package avaiable in Pip, can do gaus fitig and local maxima fiting:)
#from peakutils.plot import plot as pplot
import matplotlib.pyplot as plt
from pylab import *
import os
import sys
version = '1.4'


# In[62]:


#Sets the cell width to certan procentage
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>"))


# # Sample Selection and Read in of the Folder and Fitting Data

# In[63]:


#Folder = '26WCrY_Sa11'   # Folder/data -> where all data files which shall be used are found (and no other files)
#Folder = '27WCrY_Sa04'

#Folder = '16WCrY_Sa'
#Folder = '33WCrY'

######################  Folder = 'YReihe'
#Folder = 'NonOx'
#Folder = '34WCrY_Sa06'
#Folder = '01WCrY_3hMa'

###################### PHD Thesis
#Folder = 'W_powder'
#Folder = 'Powder_overview' # put in finished!! check double put in
#Folder = 'Bulk_Powder'     #finished

#Folder = 'LowpO2'          #finished
#Folder = 'WCr_WCrY'        #finished
#Folder = 'CrReihe'         #finished
#Folder = 'YReihe'          #finished

#Folder = 'Phasen'          #finished
#Folder = '1200C'           #finished



#Folder = 'LowPO2_bulk'     #finisched
Folder = 'Bulk_rough_smoth' #finished
#Folder = 'Bulk_thin_comp'  #finished
#Folder = 'Bulk'

CurrentDirectory = os.curdir # defines the current file where this code is saved
meas = os.listdir('%s/%s/data'%(CurrentDirectory,Folder)) # define list of data sets in the folder
dataSets = []
meas.sort() # sorting alphabetecly the file names
sec = []

for i in range (0,len(meas)): # create data arrays
    data = np.genfromtxt('%s/data/%s'%(Folder,meas[i]), skip_header = 1)
    sec.append(data[:,1].max()) ##Defining secons from the highest peak so peaks get highest peaks without overlapping
    dataSets.append(data)

print meas


# ## Import function for json database files from materialsproject.org

# In[64]:


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


# ### Test if import wokrs
# Only works if files are available on your hardrive, download from e.g. https://www.materialsproject.org/materials/mp-19005/

# In[65]:


Y6WO12 = json_data_read_reduce('daten_bank/Y6WO12_mp-19005_xrd_Cu.json')
print 'Y6WO12', len(Y6WO12), Y6WO12

WO2 = json_data_read_reduce('daten_bank/WO2_tetra_mp-19372_xrd_Cu.json')
print 'WO2', len(WO2), WO2


# ## Creating data base for oxide fitting

# # Peak Database

# In[66]:


#=======================================================================================================#
#    Peak data base for all metal and oxides: ICDD-Database (2014). PDF2014 PDF-2, www.icdd.com.        #
#    with the exception of Y6WO12 which is from: Persson, K. (2014).                                    #
#    Materials data on y6wo12 (sg:148) by materials project.                                            #
#=======================================================================================================#

W = [40.265,58.276,73.198,87.024]
Cr = [44.393, 64.5447, 81.724, 98.0349] #self fitted only if needed
Cr2O3 = [24.494,33.597,36.196, 39.749, 41.480, 44.194, 50.220, 54.852, 57.111, 58.397,
         63.449, 65.106, 72.944, 73.329, 76.851, 79.056,
         82.092, 84.239, 85.682, 86.539, 90.202, 95.328]#  33.597
WO3 = [23.083,23.707,24.099,25.956,26.587,28.776,33.331,33.640,34.022,35.525,41.524,44.880,45.354,47.228,48.431,49.326
       ,50.079,50.494,53.480,53.683,54.302,55.116]# 33.640  ,54.794,55.404]#,55.844,59.472,60.242,60.810,62.121,62.446]
#WO283 = [23.485, 23.181, 33.242, 48.077, 53.615, 59.388, 54.025]
W8O21 = [23.498, 26.038, 26.882, 31.316, 48.064, 36.000, 39.521]
WO279 = [23.512,28.214,37.044,48.093,56.694]

WO2 = json_data_read_reduce('daten_bank/WO2_tetra_mp-19372_xrd_Cu.json')

Cr2WO6 = [20.008,21.827,27.513,34.258,36.137,39.308,44.488,45.411,54.461,56.813,62.834,64.260,69.222,70.125,89.778]
YO = [29.410]#,34.038,48.888,58.075]
#CrWO4 = [19.154, 26.750, 27.681, 30.624, 36.373, 38.871, 41.089, 42.994, 44.301
#         ,50.256, 53.751, 54.794, 55.080, 57.207, 62.728, 63.736, 64.228, 68.653, 70.299, 70.663]#]
#Y2WO6 = json_data_read_reduce('daten_bank/Y2WO6_mp-565796_xrd_Cu.json')
Y6WO12 = json_data_read_reduce('daten_bank/Y6WO12_mp-19005_xrd_Cu.json')
#Y2WO6_2 = json_data_read_reduce('daten_bank/Y2WO6_2_mp-510132_xrd_Cu.json')


# filling vlaues into dic to call at key of oxide in the fitting loop
dic = { 'Cr$_2$O$_3$':sort(Cr2O3), 'WO$_3$': sort(WO3),
       'WO$_{2.79}$':sort(WO279),'Cr$_2$WO$_6$':sort(Cr2WO6),
       'W$_{8}$O$_{21}^*$':sort(W8O21),
       #,'CrWO$_4$':CrWO4
       'WO$_2$':WO2
       ,'Y$_6$WO$_{12} ^*$' : Y6WO12
      }#'Y$_2$O$_3$':[29.410,34.038,48.888,58.075]#WO283+W8O21+
#'Cr':[44.393,81.724,64.5447],#'W' : [40.265,58.276,73.198,87.024],


#dic for maker 'name': [marker,'traffic-light-color',marker ofset (0 to 0.2) to get red at top yellow, and green]
dic_marker =  { 'Cr$_2$O$_3$':['h','g',0.], 'WO$_3$': ['D','r',0.2]
               ,'WO$_{2.79}$':['^','r',0.2], 'WO$_2$':['v','r',0.2]
               , 'W$_{8}$O$_{21}^*$':['p','r',0.2] ,'Cr$_2$WO$_6$':['D','y',0.1]
               ,'Y$_2$O$_3$': ['.','g',0.],'CrWO$_4$':['*','y',0.2]
               ,'Y$_2$WO$_6 ^*$': ['d','y',0.1]
               ,'Y$_6$WO$_{12} ^*$': ['v','y',0.1]}

print dic_marker


# ## Defining ordering function for multi peak analyses

# In[67]:


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
    lambd = wavelength of the Kathode metal in angstrom default is Cu. If other source is used
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
    '''Clulating the peak position wich are really condsideren in the callulation above
    return list of values, Further more it gives back the list of the index
    
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
    '''Clulating the peak position wich should not be condsideren in the callulation above
    
    Returns
    -------
    List of values, Further more it gives back the list of the index'''
    a, b, c, d, e, f, g, h = 40.23, 44.41, 58.1, 64.6, 73.17, 81.78, 86.9, 98.04
    x = np.array(x)
    #r =  x[(a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)]
    y =   (a <x) & (x< b) ^ (c <x) & (x< d) ^ (e <x) & (x< f) ^ (g <x) & (x< h)
    y = np.logical_not(y)
    return y

#searchsort numpy based fitting ---> very fast and precise /prone if values overlap in the uncertanty range
def ox_index(a, b, uncertanty = 0.07, sort = False, unique = False, reduce_to_relevant = False):
    r'''Finds peaks in the uncertanty range (default [-0.07, 0.07]) that fit database entry
    givening back index peakpostion that doas fit the data base. 
    Base on the numpy.searchsorted function
    (https://docs.scipy.org/doc/numpy/reference/generated/numpy.searchsorted.html)
    
    Parameters
    ----------
    a : array of values (here peak postisions)
    b : array of values to compare to (here, litreture angels to expect)
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
        #b = np.sort(a)
    
    # optinal, reduce input arrays to unique values:
    if unique == True:
        a = np.unique(a)
        b = np.unique(a)
    
    # optional, reduces the input arrays, so that the longer array a values, 
    #fall in the range of the short array b
    # and visa versa
    if reduce_to_relevant == True:
        a = a[((b.min()-uncertanty)<a) & (a<(b.max()+uncertanty))]
        b = b[(a.min()<(b-uncertanty)) & ((b+uncertanty)<a.max())]
    
    #main function: finding indexes which are the intersection of the 
    #lower and higher uncertanty range of a in b
    x = np.searchsorted(a,b-uncertanty,side='left')
    y = np.searchsorted(a,b+uncertanty,side='left')-1
    index = x[x==y]
    return  index #retuns the index of the fitting parts of the array

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


# In[68]:


name = range(len(meas)) #produce list with length of measurement

for i in range(len(meas)):
    name[i] = meas[i].split('_',1)[1][:-12]#[1:3][0]#search for the first occurence of '_' splits it. take the second part and cut of -12 from the end
    print name[i]


# ## Default values for the different folders and list of exceptions

# In[69]:


#Default sizes
minpeak = None 
print minpeak
peakhight = 140.0 # hight of the resulting peak
minimal_count = 2 # minimal count of peaks for xrd phase for vergards concentration caculation
scale = peakhight*1.28
Smooth = 2 #devides and makes the the smoothing window smaller if chosen higher

# Both oxide and alloy fitting are active as default
oixdes = True
metal = True
#No oxide phase shell be fitted for [list of exeptions]
list_of_none_oxed = ['W','Cr','3.0$\,$h \nLow $p_{\mathrm{O}_2}$','0.0$\,$h','FAST7 initial'
                     ,'As produced','FAST7 befor Ox']
#where not metal phase is measured (shielded by the thick oxide layer)
list_of_none_metal = ['FAST7 rough','FAST6 humid 251$\,$h', u'CEIT \nW-10Cr-0.5Y \n10 h ox',u'FAST4']
# list of exceptions when a single additional usfull peak is to be analyse the alloy composition
exeption = [u'W-11.4Cr-0.6Y\n 3.5$\,$µm 2$\,$h Ox', u'W-11.4Cr-0.6Y\n 3.5 µm 2$\,$h Ox', 
            '2.0$\,$h', 'FAST7',u'FAST7 initial', u'FAST7 befor Ox',u'FAST7 3$\,$h 1273$\,$K',
            'FAST6',u'FAST6', u'FAST6 3523$\,$µm\n 44$\,$h Ox', '16$\,$h', '3$\,$h',
            '1.25$\,$h \n(stage 2)', '2.0$\,$h \n(stage 2)',
            u'W-Cr-Y, 3$\,$h', u'As produced', '0.1$\\,$h \n(stage 1)',
            'W-11.8Cr-0.3Y'
           ]#, '0.0$\,$h']#, u'W-Cr-Y, 3$\,$h']
# if the gausian second fittig methode is not wanted (less precize fitting)
non_gausian_metal = ['40$\,$h']

#Optional defining of names and change of default values from above if needed:
if Folder == 'WCr_WCrY':
    name[0] = 'W-11.8Cr-0.3Y' #16WCrY_Sa01
    name[1] = 'W-10.3Cr' #01WCr_Sa07
    minimal_count = 1         

if Folder == 'YReihe':
    name[0] = 'W-11.6Cr-0.4Y' #27WCrY_Sa08
    name[1] = 'W-11.4Cr-0.6Y' #33WCrY_Sa01
    name[2] = 'W-12.3Cr-0.7Y' #26WCrY_Sa07
    minimal_count = 0         # for all samples all alloy peaks are considered

if Folder == '1200C':
    name[0] = 'W-11.1Cr'#'W-11.1Cr'
    name[1] = 'W-13.5Cr-0.7Y'

if Folder == 'Powder_overview': # this one is the more default one use for other parameters as well..
    name[0] = 'W'
    name[7] = 'Cr'
    name[1] = '3$\,$h'
    name[2] = '16$\,$h'
    name[3] = '40$\,$h'
    name[4] = '60$\,$h-0.6Y'
    name[5] = '60$\,$h-1.0Y'
    name[6] = '60$\,$h-1.5Y'
    Xmin = 30.
    Xlim = 110. # maximal shown X value
    oixdes = False
    #Metal
    minpeak2 = 6.#
    peaknumber2 = 100.#65
    minimal_count = 1
    
if Folder == 'Bulk':
    name[0] = 'FAST3'
    name[1] = u'FAST4'
    name[2] = u'FAST6'
    name[3] = u'FAST7'
    #name[4] = u'CEIT \nW-10Cr-0.5Y \n10 h ox'
    minimal_count = 2
    minpeak2 = 4#6#0.08
    peaknumber2 = 80
    
if Folder == 'Bulk_Powder':
    name[0] = '60$\,$h-0.6Y'
    name[1] = 'FAST1 40$\,$h'
    name[2] = 'FAST2a'
    name[3] = 'FAST3'
    name[4] = 'FAST4 '
    name[5] = 'FAST5'
    name[6] = 'FAST6'
    name[7] = 'FAST7'
    oixdes = False
    minimal_count = 2
       
if Folder == 'W_powder':
    name[0] = 'W'
    hight = 70
    oixdes = False
    
if Folder == 'Phasen':
    name[0] = u'0.0$\,$h'#
    name[1] = '0.1$\,$h \n(stage 1)'
    name[2] = '1.25$\,$h \n(stage 2)'
    name[3] = '2.0$\,$h \n(stage 2)'
    name[4] = '3.0$\,$h \n(stage 3)'
    name[5] = '8.0$\,$h \n(stage 3)'
    minimal_count = 1

if Folder == 'LowpO2':
    name[2] = u'As produced'
    name[1] = u'W-Cr, 3$\,$h'
    name[0] = u'W-Cr-Y, 3$\,$h'
    oixdes = False
    minimal_count = 2

if Folder == 'CrReihe':
    name[0] = 'W-8.0Cr-0.6Y'  #30WCrY_Sa09
    name[1] = 'W-10.2Cr-0.6Y' #32WCrY_Sa04
    name[2] = 'W-11.4Cr-0.6Y' #33WCrY_Sa04
    name[3] = 'W-13.5Cr-0.7Y' #34WCrY_Sa03
    minimal_count = 0 #every peak is considered, here one have be carefull in the peak detection setting
    
if Folder == 'LowPO2_bulk':
    name[0] = u'FAST7 initial'
    name[1] = u'FAST7 3$\,$h 1273$\,$K'
    minpeak2 = 4#6#0.08
    peaknumber2 = 80
    oixdes = False
    minimal_count = 1

if Folder == 'Bulk_thin_comp':
    name[2] = u'W-11.4Cr-0.6Y\n 3.5$\,$µm 2$\,$h Ox'
    name[1] = u'W-13.5Cr-1.1Y\n 7.5$\,$µm 6$\,$h Ox'
    name[0] = u'FAST6 3523$\,$µm\n 44$\,$h Ox'
    minimal_count = 2
    minpeak2 = 4#6#0.08
    peaknumber2 = 80
    
if Folder == 'Bulk_rough_smoth':
    name[0] = u'FAST7 befor Ox'
    name[1] = 'FAST7 standard'# 44 h Ox'
    name[2] = 'FAST7 rough'# 44 h Ox'
    #oxide
    minpeak = 0.045#0.08
    peaknumber = 200
    minpeak = [0.03,0.02,0.06]#*len(name)
    peaknumber = [150]*len(name)
    #metal
    minpeak2 = [10.]*len(name)#8.3#4#8.1##0.06#0.08 #percent of maximal peak
    peaknumber2 = [100.]*len(name)#65

if oixdes == True:
    start = 20
    stop = 100
else:
    start = 39
    stop = 95
print 'SETTINGS:\n'
print Folder, name
print 'start-stop',start, stop


# ## Plot Settings

# In[70]:


get_ipython().magic(u"config InlineBackend.figure_format='svg'")
# Update the matplotlib configuration parameters:
matplotlib.rcParams.update({'font.size': 14, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
matplotlib.rcParams['svg.fonttype'] = 'none'

from matplotlib import rc
rc('lines', lw=1.5)
if len(name)>8:
    color= cm.tab10(np.linspace(0,1,len(name)))#Dark2
else:
    color= cm.Dark2(np.linspace(0,1,len(name)))
    #color= cm.inferno(np.linspace(0,1,len(name)+2))
    #color = color[::-1]#reversed for comparison to other oxidation graph
    
#Automatic figure hight setting by number of lines
print len(name)
if len(name)<5:
    hight = len(name)*2.
elif len(name)<7:
    hight = len(name)*1.35
else:
    hight = len(name)*1

rc('figure', figsize=(6.2, hight))

print hight


# # Peak fitting
# ## Calculate Cr content in W-Cr phase from peakshift relative to the tungsten peak

# In[71]:


#Default values
aCr = 2.8846#2.8846#2.885 # lattice constant chromium in angstrom 288.46 
aW = 3.16  # lattice constant tungsten in angstrom
lambd = 1.54184 # wavelength of the x-ry Cu source in angstrom
atomicMassW = 183.84
atomicMassCr = 51.9961


# ## Check if fitt settings are aproriet for the sample investigated
# ### Change mindpeak and peaknumber variable for oxide fitting
# For different peak fitting parameter, list can be filled
# Current folder and samples:

# In[72]:


print Folder
print name


# # Oxides

# In[73]:


#default
minpeak = [0.05]*len(meas)
peaknumber = [150]*len(meas)
#HERE set your wanted values if default are not sufficient

if Folder == 'Bulk_rough_smoth':
    minpeak = [0.02,0.03,0.05]#*len(meas)
if Folder == 'Bulk_thin_comp':
    minpeak = [0.03,0.05,0.04]#*len(meas)
if Folder == 'Bulk':
    minpeak = [0.03]*len(meas)
    peaknumber = [150]*len(meas)
if Folder == 'Phasen':
    minpeak = [0.05]*len(meas)
    minpeak[3] = 0.08
if Folder == 'CrReihe':
    minpeak[3], minpeak[2], minpeak[0] = 0.02, 0.03, 0.04
if Folder == 'WCr_WCrY':
    minpeak[0] = 0.03
 
Smooth = 2
if oixdes == True:
    for i in reversed(range(0,len(meas))):
        x_all = dataSets[i][:,0]
        x = x_all[(start < x_all) & (x_all < stop)]
    
    
        y = peakhight/sec[i]*(dataSets[i][(start <x_all) & (x_all<stop),1])+scale*i
        ## smoothing the data for the fit, makes fiitig more reliable:
        window = (len(x)/(stop-start))/Smooth
        window = window if window%2 != 0 else window+1#checks is window is odd
        y_smooth = savgol_filter(y, window_length=window, polyorder=3)
        print name[i],len(x)/(stop-start),u'pixel/°', 'window =',window, minpeak[i]
    
        plot(x,y_smooth,c=color[i],lw=1,label='smoothed data')
        plot(x,y,c=color[i],lw=0.5,color='red')
    
        scatter(x,y,s=2,c='red',label='data' if 'data' not in plt.gca().get_legend_handles_labels()[1] else '')
        #plt.axvline(x=26.16878048)
    
        fill(x, y, zorder=10,color=color[i],alpha=0.2)
        plt.text(stop,peakhight/sec[i]*(dataSets[i][100,1])+scale*i, name[i]+' [{}] {:.2f}'.format(i,minpeak[i]), color=color[i])#prints the name to the line
        #fitiing and ploting
        indexes = peakutils.indexes(y, thres=minpeak[i], min_dist=len(x)/peaknumber[i]) #finding peak index position /max(y)
        print type(indexes),indexes, '\n', x[indexes]
        peaks_x = peakutils.interpolate(x, y, ind=indexes, width = window)#gausian fit better presition
        peaks_x = np.where(np.abs(peaks_x-x[indexes])>1.,x[indexes],peaks_x)#check wrong fit
        #print name[i],peaks_x[6], y[ox_index(y[indexes],y)]
        
    

        #check if at given x value y value fits the data
    
        #x[Ox_index(peaks_x,x,uncertanty=[-0.02,0.02])]
        #scatter(x[Ox_index(peaks_x,x,uncertanty=[-0.02,0.02])],y[Ox_index(peaks_x,x,uncertanty=[-0.02,0.02])])
        if (oixdes==True):# and (name[i] not in list_of_none_oxed)
            scatter(peaks_x,y[indexes]+0.15*peakhight)
        #print name[i], peaks_x

    xlim(start, stop)
print Folder


# # Metal/Alloy
# ### Change mindpeak2 and peaknumber2 variable for alloy phase fitting

# In[74]:


minpeak2 = [12]*len(name)#8.3#4#8.1##0.06#0.08 #percent of maximal peak
peaknumber2 = [90.]*len(name)#65

if Folder == 'Bulk_rough_smoth':
    minpeak2 = [2.1,2,2]#*len(name)
    peaknumber2 = [80.]*len(name)#65
if Folder == 'LowPO2_bulk':
    minpeak2 = [2, 2]#*len(name)# [1] 1.0 #percent of maximal peak
    peaknumber2 = [120.]*len(name)#65
if Folder == 'Bulk_thin_comp':
    minpeak2 = [8,14,10]#*len(name)#8.3#4#8.1##0.06#0.08 #percent of maximal peak
    peaknumber2 = [90.,110.,90.]#*len(name)#65
if Folder == 'Phasen':
    minpeak2[5],minpeak2[4], minpeak2[3], minpeak2[1], minpeak2[0] = 5, 3.5, 13, 7.0, 16
if Folder == 'Bulk_Powder':
    minpeak2 = [4.4]*len(name)#8.3#4#8.1##0.06#0.08 #percent of maximal peak
    minpeak2[0]= 5
    Smooth = 2 #deider for window for smoothing
if Folder == 'Powder_overview':
    minpeak2 = [6]*len(name)#8.3#4#8.1##0.06#0.0,8 #percent of maximal peak
    minpeak2[1], minpeak2[2], minpeak2[3], minpeak2[6] , minpeak2[7] = 2, 3.5, 5, 8,4
    peaknumber2[3] = 80
if Folder == 'Bulk':
    minpeak2 = [6]*len(meas)
if Folder == 'LowpO2':
    minpeak2 = [6.5]*len(name)#8.3#4#8.1##0.06#0.08 #percent of maximal peak
    minpeak2[2] = 16
    #Smooth = 3 #??
if Folder == 'CrReihe':
    minpeak2 = [12]*len(name)
if Folder == 'YReihe':
    minpeak2[0] = 6
if Folder == 'WCr_WCrY':
    minpeak2[1] = 6
    minpeak2[0] = 4
    
print minpeak2

#color = color[bla]
for i in reversed(range(0,len(meas))):
    x_all = dataSets[i][:,0]
    x = x_all[(start < x_all) & (x_all < stop)]
    y = peakhight/sec[i]*(dataSets[i][(start <x_all) & (x_all<stop),1])+scale*i
    
    ##smoothing data
    window = (len(x)/(stop-start))/Smooth
    window = window if window%2 != 0 else window+1 # checks is window is odd
    print 'window',window
    y_smooth = savgol_filter(y, window_length=window, polyorder=3)
    plot(x,y,c=color[i])
    fill(x, y, zorder=10,color=color[i],alpha=0.2)
    plt.text(stop,peakhight/sec[i]*(dataSets[i][100,1])+scale*i, name[i]+'\t[{}] {:.1f}'.format(i,minpeak2[i]), color=color[i])
    #prints the name to the line
    
    
    if (metal == True):# and (name[i] not in list_of_none_metal):#can be turned on if filter list should apply
        indexes = peakutils.indexes(y_smooth, thres=minpeak2[i]/peakhight, min_dist=len(x)/peaknumber2[i])
        #finding peak index position /max(y)
        peaks_x2 = peakutils.interpolate(x, y_smooth, ind=indexes, width=window)#gausian refit for better precision
        if name[i] in non_gausian_metal:
             peaks_x2 = x[indexes]
        peaks_x2 = np.where(np.abs(peaks_x2-x[indexes])>1.,x[indexes],peaks_x2)#check wrong fit
        aWCr = peak(peaks_x2)#calulating the lattice constant for the bcc phase
        x_reduced, yin = relevant(peaks_x2)# reducing all peaks to the relevant bcc ones  (here just metal)
        y_reduced = y_smooth[indexes[yin]]#same for peakhight
        colorX = color[i]
        scatter(x_reduced,y[indexes[yin]]+0.15*peakhight, marker= '+')
        #scatter(x_reduced,y[indexes[yin]]+0.15*peakhight)
        
        #print name[i], aWCr # print name and lattice constant for given peak
        
xlim(start,stop)
#xlim(60,90)
#xlim(70, 80)#, stop)


# # The Plot

# In[75]:


### dic_peaks = {} #open new dict where the name and peak data is filled in
fig1 = plt.figure(facecolor='white')
ax1 = plt.axes(frameon=False)
ax1.set_frame_on(False)
ax1.get_xaxis().tick_bottom()
ax1.axes.get_yaxis().set_visible(False)#True for visible Y achses
ax1.axes.get_xaxis().set_visible(True)
xlabel(u'$2\,\Theta$ [°]')
#ylabel(u'Intensity [a.u.]')

#import Plot_XRD #imports the look of the plot no boder and the right x and y labels...
dic_peaks = {}

for i in reversed(range(0,len(meas))):
    x_all = dataSets[i][:,0]
    x = x_all[(start < x_all) & (x_all < stop)]
    y=peakhight/sec[i]*(dataSets[i][(start <x_all) & (x_all<stop),1])+scale*i
    ###################producing the lines area for The W-Cr BCC systems#################
    if (metal== True) & (oixdes==True):
        comun_shift = 0
    else:
        comun_shift = 1
    if metal== True:
        if i == 0:
            for j in range(len(W)):
                plt.vlines(W[j],0,(scale)*len(name)-6, lw=0.5,alpha=0.7)
                plt.axvspan(W[j], Cr[j] if Cr[j]<stop else stop, alpha=0.0125, color='black')
            plt.text(W[-1-comun_shift],(scale)*len(name),'W',fontsize = 12,horizontalalignment='center')
            for j in range(len(Cr)-comun_shift):
                plt.vlines(Cr[j],0,(scale)*len(name)-6, lw=0.5,alpha=0.4)
            plt.text(Cr[-1-comun_shift],(scale)*len(name),'Cr',fontsize = 12,horizontalalignment='center')
            plt.text(W[-1-comun_shift]+(Cr[-1-comun_shift]-W[-1-comun_shift])/2,(scale)*len(name),                     u'\u21C4',fontsize = 13,horizontalalignment='center')
        
    plot(x,y,c=color[i])
    fill(x, y, zorder=10,color=color[i],alpha=0.2)
    plt.text(stop,peakhight/sec[i]*(dataSets[i][100,1])+scale*i, name[i], color=color[i])#prints the name to the line
    print name[i]
    ## smoothing the data for the fit, makes fiitig more reliable:
    window = (len(x)/(stop-start))/Smooth   
    window = window if window%2 != 0 else window+1#checks is window is odd
    y_smooth = savgol_filter(y, window_length=window, polyorder=3)
    
    #==================================================================================#
    #              FIT/Plot concentrations... Vegard approximation                     #
    #==================================================================================#
    if (metal == True) and (name[i] not in list_of_none_metal):
        ##making the fit on the smoothed data
        indexes = peakutils.indexes(y_smooth, thres=minpeak2[i]/peakhight, min_dist=len(x)/peaknumber2[i])
        #finding peak index position /max(y)
        peaks_x2 = peakutils.interpolate(x, y_smooth, ind=indexes, width=window)#gausian refit for better precision
        if name[i] in non_gausian_metal:
             peaks_x2 = x[indexes]
        peaks_x2 = np.where(np.abs(peaks_x2-x[indexes])>1.,x[indexes],peaks_x2)#check wrong fit
        aWCr = peak(peaks_x2)#calulating the lattice constant for the bcc phase
        x_reduced, yin = relevant(peaks_x2)# reducing all peaks to the relevant bcc ones  (here just metal)
        y_reduced = y[indexes[yin]]#same for peakhight
        colorX = color[i]
        inds =[] #pruce empty list to fill with indices for each bin (histogram) idividualy
        print len(aWCr),'{:.2}'.format(aWCr.std()),'standart deviation lattice constant'
        
        if aWCr.std()>0.004: #filtering if there are more than one cubic phase, only shown if 3 or more## turn off by increasing standdeviation
            #check length of a, and find proper spacing of bins
            print 'ON!!length{}'.format(len(aWCr)), sort(aWCr)#tells u that filtering is active
            #==================================================================================#
            #    increasing the number of bins till a reasonable amount of phases is found     #
            #==================================================================================#
            print 'frequency of phase peak',plane_count(peaks_x2)[0]
            print 'minimal number of phases',plane_count(peaks_x2)[1]
            bin_num = plane_count(peaks_x2)[1]
            while True:
                #creats bins which are used to filter the data (equivalent to a histogram)
                bins = np.linspace(aWCr.min(),aWCr.max()+0.00000001,bin_num)
                print 'len(bins)',len(bins), bins
                #creats temporary emtpy np.array in the right size
                counts = np.zeros_like(bins)
                # find the appropriate bin foreach a, right means right vlaue<bin index>bin-index
                k = np.searchsorted(bins, aWCr,side='right')
                # add one for every fitting value to a bin
                np.add.at(counts, k, 1)
                #counts maximal count
                frequency = counts.max()
                print 'frequency', frequency
                print counts
                # breaks if condition 1 and 2 is fullfilled
                if (frequency <= plane_count(peaks_x2)[0]) & ((counts>0).sum()  >= plane_count(peaks_x2)[1]):
                    break
                #increase bin number for better separetion
                bin_num += 1
            
            #creats array of lenght bins with the aproperiat index values of bins
            v_ = np.arange(len(bins))
            
            #checks at which bin index the count is higher than minmal_count 2(default or 1 for thinfilm samples) (so significant)
            vin =  v_[(counts>2)]
            
            #check if any True value is found
            if (((counts>2).sum() == False) or (minimal_count < 2)):
                vin =  v_[(counts>(1 if minimal_count == 2 else minimal_count))]#if not, reduce condition to more than 1
                print 'FILTERT with>{}'.format(1 if minimal_count == 2 else minimal_count), (counts>(1 if minimal_count == 2 else minimal_count))
            if name[i] in exeption:
                vin =  v_[(counts>0)]#if not, reduce condition to more than 1
                print 'FILTERT expection with=1,', (counts>0),vin
            for n in range(len(vin)):
                inds.append([k==vin[n]])    
        else:
            inds.append([aWCr==aWCr])
            vin  = list([True])
            print 'only one phase'
        #print 'inds',inds
        for n in reversed(range(len(vin))):
            atPerCr = 1.-(aWCr[inds[n]]-aCr)/(aW-aCr)#calculate at% Cr
            wtPerCr = atPerCr*atomicMassCr/(atomicMassW*(1-atPerCr)+atPerCr*atomicMassCr)# wt%Cr
            scatter(x_reduced[inds[n]],y_reduced[inds[n]]+0.18*peakhight, s=40, marker='s'                    ,color = color[i], alpha = len(vin)/(len(vin)+float(n))                    ,label = '{}$\cdot$W-{:.0f}Cr'.format(len(wtPerCr),abs(wtPerCr.mean())*100))

    #===============================================================================================#
    #                      Fitting data from libary (pure Cr,W and oxides)                          #
    #===============================================================================================#
    indexes = peakutils.indexes(y_smooth, thres=minpeak[i], min_dist=len(x)/peaknumber[i]) # , fontsize = 16 finding peak index position /max(y)
    peaks_x = peakutils.interpolate(x, y_smooth, ind=indexes, width=window)#gausian fit better presition
    #check if any peaks were lost and if get rid of them
    peaks_x = np.where(np.abs(peaks_x-x[indexes])>0.8,x[indexes],peaks_x)#checks if fit is off after interpolation by a sustancial amount.
    dic_peaks[name[i]] = peaks_x
    rel = peaks_x==peaks_x#takes all, in the list into consideration
    #recues double fitting in the alloy zones, can be turned of for specific samples by adding to the list_of_none_metal
    if (metal == True) & (name[i] not in list_of_none_metal):
        rel = irrelevant(peaks_x)
        peaks_x = peaks_x[rel]
     
    # ad a minmal diff function to the original so no double fitting...?? possible?
    if (oixdes==True) and (name[i] not in list_of_none_oxed):
        #print name[i], (name[i] not in list_of_none_oxed)
        for key in dic:
            temp_index_ox = ox_index(peaks_x,dic[key])#, uncertanty=0.12) #<---- if comparison presicion should be reduced (not higher than 0.15)
            if Folder == 'Bulk':#exception for  bulk
                temp_index_ox = ox_index(peaks_x,dic[key], uncertanty=0.15)
            #print key,len(peaks_x),
            dic_array = np.array(dic[key])
            peak_array = peaks_x[temp_index_ox]
            print key, len(temp_index_ox), '{:.3f}'.format(np.abs(peak_array[ox_index(peak_array,dic[key])]-dic_array[ox_index(dic[key],peak_array)]).mean())
            if (len(temp_index_ox) >= 2):#<----------checks if array is empty=1, or filters for more than one peak if set to 2
                #checks the average distance fit position and databank
                if (np.abs(peak_array[ox_index(peak_array,dic_array)]-dic_array[ox_index(dic_array,peak_array)]).mean()<0.056):
                    scatter(peaks_x[temp_index_ox],y[indexes][rel][temp_index_ox]+(0.1+dic_marker[key][2])*peakhight                            ,marker=dic_marker[key][0],color=dic_marker[key][1],s = 40,alpha=0.8                            ,label = key if key not in plt.gca().get_legend_handles_labels()[1] else '')
                             #label = checks if the key is already in the list of labels and so avoids double naming of labels
    print '\n'
dic_oxtest = dic_peaks

plt.legend(handlelength=0, bbox_to_anchor=(.05, 0.98), ncol=1, mode=None)#left

print Folder
#For last plot make '_' for testing '_test'
test = '_test'
print '%s/06_XRD%snew%s'%(Folder,test,Folder)

savefig('%s/06_XRD%snew%s.svg'%(Folder,test,Folder), bbox_inches = 'tight', transparent=True)#saves the fig in chooesen format
savefig('%s/06_XRD%snew%s.pdf'%(Folder,test,Folder), bbox_inches = 'tight', transparent=True)#bbox_inches ='
savefig('%s/06_XRD%snew%s.png'%(Folder,test,Folder), bbox_inches = 'tight', transparent=True,dpi = 300)#bbox_inches ='


# # End
# Here the fitting parameters can be tested if they are reasonably well for the Folder

# In[59]:


plane_count(peaks_x2)
print 'maximal_frequency', plane_count(peaks_x2)[0]
print 'minimalbin', plane_count(peaks_x2)[1]


# In[71]:


#rouding based Fitting---> slower but little more stable
def index_round(peaks_x,ox,dec = 1):
    peaks_x = np.array(peaks_x)
    ox = np.array(ox)
    a = np.round(peaks_x, decimals = dec)
    b = np.round(ox, decimals = dec)
    c = np.isin(a,b)
    d = np.isin(b,a)
    #if len(c) > len(d):
     #   c = np.searchsorted(a,a[c])
    #else:
    #    c = np.searchsorted(a[c],b[d],side='left')
    return c


# ## reading in data directly from a database (not finished)

# In[72]:


import requests

def get_api_key(api_key):
    '''Returns the string stored in the file with the path "api_key" '''
    try:
        with open(api_key) as inputFileHandle:
            return inputFileHandle.read()[:-1]#-1 to get ride of the new line '\n'
    except IOError:
        sys.stderr.write('api_key - Error: Could not open {}, check if it is there?\n'.format(api_key))
        sys.exit(-1)


# In[73]:


r = requests.get('https://www.materialsproject.org/rest/v1/materials/Fe2O3/vasp', params= {'API_KEY': get_api_key('daten_bank/key')})


# In[ ]:


#r.content

