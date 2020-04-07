
from scipy import *
from numpy import *
from math import *
from numpy import matrix

import copy
import matplotlib.pyplot as pl
import numpy as np
import matplotlib.backends.backend_pdf
from matplotlib import pyplot
#import pylab

import sys
from sys import exit
from math import pi
import calibfunctions as cf

folder = str(sys.argv[1])

if "-h" in folder:
    print "Usage of fastAnalysis.py:"
    print "Uses Python version 2.7. check currently loaded python with module list or module avail"
    print "command format: \n\n python fastAnalysis.py {folder} {extension} \n"
    print "{folder} = folder to search for the Laser output files of the format 'tenkLaserRun_[Position]-[wavelength]-{extension}.txt'"
    
    quit

incode = str(sys.argv[2])

inplace = [""]*6
inplace[0] = "Top-532"
inplace[1] = "Middle-532"
inplace[2] = "Bottom-532"
inplace[3] = "Top-213"
inplace[4] = "Middle-213"
inplace[5] = "Bottom-213"

times = [[]]*6
cutms = [[]]*6

rowratios = []
ucratios = []
ucstring=incode+"\n"

for n in range(0,6):    
    rrats = []
    ucrat = 0
    times[n],cutms[n] = cf.fastprocessfile(folder+"/tenkLaserRun10x_"+inplace[n]+"-"+incode+".txt")
    rrats = cf.rowRatios(times[n],cutms[n],(inplace[n]+"-offoff"+incode),2)
    ucrat = cf.wecalcRatios(times[n],cutms[n])

    ucratios.append(ucrat)
    ucstring = ucstring+inplace[n]+" u/c= "+str(ucrat)+" \n"
    for i in range(0,len(rrats)):
        rowratios.append(rrats[i])

simfig,sqrfig,nabfig = cf.figureOfMerit(rowratios,ucratios,True)
ucstring = ucstring+" simple = "+str(simfig)+"\n square = "+str(sqrfig)+"\n nonabs = "+str(nabfig)+"\n\n\n"


f = open("fastRatios.txt","a")
f.write(ucstring)
f.close()
    
