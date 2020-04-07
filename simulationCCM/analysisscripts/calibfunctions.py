from scipy import *
from numpy import *
from math import *
from numpy import matrix

import copy
import matplotlib.pyplot as pl
import numpy as np
import matplotlib.backends.backend_pdf
from matplotlib import pyplot
from scipy.stats import norm
#import pylab                                                                                              

import sys
from sys import exit
from math import pi

#This function build a plot of a DAQ window for a given event. Entries are double[] data=PE times, string[] coated=PMT name for each of those PE times, string title=plot title, string yval=plot y label, and int shift=how to move the 0 time of the event (+100=event starts 100 ns late).
def buildwindow(data,pmt,title,yval,shift):
    #these values can be changed to limit the DAQ window shown (currently: .95 to 1.35 microseconds after '0' value
    xmin=950
    xmax=1350
    nbins=np.linspace(xmin,xmax,201)

    #Make the shift random or standardized
    #    shift = np.random.uniform(1000,5000)
    #    shift = 1000 
    
    #establish some variables
    nevents = len(data)
    data2=[]
    data3=[]
    ws=[]
    ws3=[]
    data2 = data[:]
    ws = [1./30.]*len(data2)

    #converts each PE into a triangular waveform 20 ns wide, spreading the PE value across all of it
    for i in range(0,nevents):
        data2[i] = data2[i]+shift
        data2.append(data2[i]+18)
        ws.append(float(1./30.))
        for n in range(1,5):
            nv = n+0.001-0.001
            data2.append(data2[i]+2*n)
            ws.append(float(n+1)/30.)
            data2.append(data2[i]+(18-2*n))
            ws.append(float(n+1)/30.)
            
        #Creates a seperate waveform for the uncoated data
        if "coated" not in pmt[i]:
            data3.append(data2[i])
            data3.append(data2[i]+18)
            ws3.append(float(1./30.))
            ws3.append(float(1./30.))
            for n in range(1,5):
                data3.append(data2[i]+2*n)
                ws3.append(float(n+1)/30.)
                data3.append(data2[i]+(18-2*n))
                ws3.append(float(n+1)/30.)

    #create a histogram of the PE values per each bin (2ns wide) as in the actual data
    xhist,xedges = histogram(data2,bins=nbins,weights=ws,range=((xmin),(xmax)))
    xhist2,xedges2 = histogram(data3,bins=nbins,weights=ws3,range=((xmin),(xmax)))

    #plot that waveform. Will be saved to a file if done outside this function
    x = array( [ (xedges[i]+xedges[i+1])/2. for i in range(len(xedges)-1) ] )
    x2 = array( [ (xedges2[i]+xedges2[i+1])/2. for i in range(len(xedges2)-1) ] )
    pl.errorbar(x,xhist,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x2,xhist2,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.ylabel(yval,fontsize=20)
    pl.title(title,fontsize=24)

#Same as above, but has a user input window (from 0 to xmax).
def buildfullwindow(data,pmt,title,yval,shift,xmax):
    xmin=0
    xcount = int(xmax/2+1)
    nbins=np.linspace(xmin,xmax,xcount)
    #    shift = np.random.uniform(1000,5000)
    #    shift = 1000 
    nevents = len(data)
    data2=[]
    data3=[]
    ws=[]
    ws3=[]
    data2 = data[:]
    
    ws = [1./30.]*len(data2)
    for i in range(0,nevents):
        data2[i] = data2[i]+shift
        data2.append(data2[i]+18)
        ws.append(float(1./30.))
        for n in range(1,5):
            nv = n+0.001-0.001
            data2.append(data2[i]+2*n)
            ws.append(float(n+1)/30.)
            data2.append(data2[i]+(18-2*n))
            ws.append(float(n+1)/30.)
            
        if "coated" not in pmt[i]:
            data3.append(data2[i])
            data3.append(data2[i]+18)
            ws3.append(float(1./30.))
            ws3.append(float(1./30.))
            for n in range(1,5):
                data3.append(data2[i]+2*n)
                ws3.append(float(n+1)/30.)
                data3.append(data2[i]+(18-2*n))
                ws3.append(float(n+1)/30.)

    xhist,xedges = histogram(data2,bins=nbins,weights=ws,range=((xmin),(xmax)))
    xhist2,xedges2 = histogram(data3,bins=nbins,weights=ws3,range=((xmin),(xmax)))

    x = array( [ (xedges[i]+xedges[i+1])/2. for i in range(len(xedges)-1) ] )
    x2 = array( [ (xedges2[i]+xedges2[i+1])/2. for i in range(len(xedges2)-1) ] )
    pl.errorbar(x,xhist,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x2,xhist2,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.ylabel(yval,fontsize=20)
    pl.title(title,fontsize=24)

#Builds an accumulated DAQ waveform. Also does so for multiple sets of data (3 currently). Most recent use was building the three waveforms for the three laser positions. Also saves their relative ratios to a text file. Adds the fig and fn inputs (as well as multiple copies of data and pmt), which are respectively figure fig=the figure to be added to and int fn=the file to be saved to, see below for details
def buildaccumulated(data,pmt,datam,pmtm,datab,pmtb,title,yval,fig,fn):
    ax=fig.add_subplot(1,1,1)
    
    xmin=950
    xmax=1150
    nbins=np.linspace(xmin,xmax,101)
#    shift = np.random.uniform(1000,5000)
    shift = 1000 
    nevents = len(data)
    data2=[]
    data3=[]
    ws=[]
    ws3=[]
    data2 = data[:]

    ws = [1./30.]*len(data2)
    for i in range(0,nevents):
        data2[i] = data2[i]+shift
        data2.append(data2[i]+18)
        ws.append(float(1./30.))
        for n in range(1,5):
            nv = n+0.001-0.001
            data2.append(data2[i]+2*n)
            ws.append(float(n+1)/30.)
            data2.append(data2[i]+(18-2*n))
            ws.append(float(n+1)/30.)

        if "coated" not in pmt[i]:
            data3.append(data2[i])
            data3.append(data2[i]+18)
            ws3.append(float(1./30.))
            ws3.append(float(1./30.))
            for n in range(1,5):
                data3.append(data2[i]+2*n)
                ws3.append(float(n+1)/30.)
                data3.append(data2[i]+(18-2*n))
                ws3.append(float(n+1)/30.)

    neventsm = len(datam)
    datam2=[]
    datam3=[]
    wsm=[]
    wsm3=[]
    datam2 = datam[:]

    wsm = [1./30.]*len(datam2)
    for i in range(0,neventsm):
        datam2[i] = datam2[i]+shift
        datam2.append(datam2[i]+18)
        wsm.append(float(1./30.))
        for n in range(1,5):
            nv = n+0.001-0.001
            datam2.append(datam2[i]+2*n)
            wsm.append(float(n+1)/30.)
            datam2.append(datam2[i]+(18-2*n))
            wsm.append(float(n+1)/30.)

        if "coated" not in pmtm[i]:
            datam3.append(datam2[i])
            datam3.append(datam2[i]+18)
            wsm3.append(float(1./30.))
            wsm3.append(float(1./30.))
            for n in range(1,5):
                datam3.append(datam2[i]+2*n)
                wsm3.append(float(n+1)/30.)
                datam3.append(datam2[i]+(18-2*n))
                wsm3.append(float(n+1)/30.)

    neventsb = len(datab)
    datab2=[]
    datab3=[]
    wsb=[]
    wsb3=[]
    datab2 = datab[:]

    wsb = [1./30.]*len(datab2)
    for i in range(0,neventsb):
        datab2[i] = datab2[i]+shift
        datab2.append(datab2[i]+18)
        wsb.append(float(1./30.))
        for n in range(1,5):
            nv = n+0.001-0.001
            datab2.append(datab2[i]+2*n)
            wsb.append(float(n+1)/30.)
            datab2.append(datab2[i]+(18-2*n))
            wsb.append(float(n+1)/30.)

        if "coated" not in pmtb[i]:
            datab3.append(datab2[i])
            datab3.append(datab2[i]+18)
            wsb3.append(float(1./30.))
            wsb3.append(float(1./30.))
            for n in range(1,5):
                datab3.append(datab2[i]+2*n)
                wsb3.append(float(n+1)/30.)
                datab3.append(datab2[i]+(18-2*n))
                wsb3.append(float(n+1)/30.)
                
    xhist,xedges = histogram(data2,bins=nbins,weights=ws,range=((xmin),(xmax)))
    xhist2,xedges2 = histogram(data3,bins=nbins,weights=ws3,range=((xmin),(xmax)))
    xhistm,xedgesm = histogram(datam2,bins=nbins,weights=wsm,range=((xmin),(xmax)))
    xhistm2,xedgesm2 = histogram(datam3,bins=nbins,weights=wsm3,range=((xmin),(xmax)))
    xhistb,xedgesb = histogram(datab2,bins=nbins,weights=wsb,range=((xmin),(xmax)))
    xhistb2,xedgesb2 = histogram(datab3,bins=nbins,weights=wsb3,range=((xmin),(xmax)))
    
    x = array( [ (xedges[i]+xedges[i+1])/2. for i in range(len(xedges)-1) ] )
    x2 = array( [ (xedges2[i]+xedges2[i+1])/2. for i in range(len(xedges2)-1) ] )
    pl.errorbar(x,xhist,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x2,xhist2,yerr=0,drawstyle='steps-mid',capsize=0)
    xm = array( [ (xedgesm[i]+xedgesm[i+1])/2. for i in range(len(xedgesm)-1) ] )
    xm2 = array( [ (xedgesm2[i]+xedgesm2[i+1])/2. for i in range(len(xedgesm2)-1) ] )
    pl.errorbar(xm,xhistm,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(xm2,xhistm2,yerr=0,drawstyle='steps-mid',capsize=0)
    xb = array( [ (xedgesb[i]+xedgesb[i+1])/2. for i in range(len(xedgesb)-1) ] )
    xb2 = array( [ (xedgesb2[i]+xedgesb2[i+1])/2. for i in range(len(xedgesb2)-1) ] )
    pl.errorbar(xb,xhistb,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(xb2,xhistb2,yerr=0,drawstyle='steps-mid',capsize=0)
    
    nett = sum(xhist2)
    netm = sum(xhistm2)
    netb = sum(xhistb2)
    pl.ylabel(yval,fontsize=20)
    pl.title(title,fontsize=24)

    textstr = '$N_t$=%5.0f \n $N_m$=%5.0f \n $N_b$=%5.0f \n $Nu_t$=%3.3f \n $Nu_m$=%3.3f \n $Nu_b$=%3.3f'% (nevents,neventsm,neventsb,nett,netm,netb)
    props=dict(boxstyle='round',facecolor='white',alpha=0.5)
    pl.text(0.75,0.95,textstr,fontsize=16,verticalalignment='top',transform=ax.transAxes,bbox=props)
    
    #potential filenames for saving. Add more if you want.
    filename = ["LaserRatios.txt","rerunLaserRatios.txt"]
    f = open(filename[fn],"a")
    f.write(title)
    f.write('\n B/T = %3.3f \n M/T = %3.3f \n uB/uT = %3.3f \n uM/uT = %3.3f \n uB/cB = %3.3f \n uM/cM =  %3.3f \n uT/cT = %3.3f \n'% ((float(neventsb)/nevents),(float(neventsm)/nevents),(float(netb)/nett),(float(netm)/nett),(float(netb)/24/(float(neventsb-netb)/96)),(float(netm)/24/(float(neventsm-netm)/96)),(float(nett)/24/(float(nevents-nett)/96))))
    f.close()
    
#function to build the row ratios. takes in the same values as above, just in a different order.
def buildrows(data,pmt,title,yval,fig,shift,fn):
    ax=fig.add_subplot(1,1,1)

    filename = ["LaserRatios.txt","rerunLaserRatios.txt"]
    f = open(filename[fn],"a")
    f.write(title)

    xmin=950
    xmax=1150
    nbins=np.linspace(xmin,xmax,101)
#    shift = np.random.uniform(1000,5000)
    nevents = len(data)
    data1=[]
    data2=[]
    data3=[]
    data4=[]
    data5=[]
    datan=[]
    ws1=[]
    ws2=[]
    ws3=[]
    ws4=[]
    ws5=[]
    datan = data[:]

    for i in range(0,nevents):
        datan[i] = datan[i] + shift
        if "R1" in pmt[i]:
            data1.append(datan[i])
            data1.append(datan[i]+18)
            ws1.append(float(1./30.))
            ws1.append(float(1./30.))
            for n in range(1,5):
                data1.append(datan[i]+2*n)
                ws1.append(float(n+1)/30.)
                data1.append(datan[i]+(18-2*n))
                ws1.append(float(n+1)/30.)
        if "R2" in pmt[i]:
            data2.append(datan[i])
            data2.append(datan[i]+18)
            ws2.append(float(1./30.))
            ws2.append(float(1./30.))
            for n in range(1,5):
                data2.append(datan[i]+2*n)
                ws2.append(float(n+1)/30.)
                data2.append(datan[i]+(18-2*n))
                ws2.append(float(n+1)/30.)
        if "R3" in pmt[i]:
            data3.append(datan[i])
            data3.append(datan[i]+18)
            ws3.append(float(1./30.))
            ws3.append(float(1./30.))
            for n in range(1,5):
                data3.append(datan[i]+2*n)
                ws3.append(float(n+1)/30.)
                data3.append(datan[i]+(18-2*n))
                ws3.append(float(n+1)/30.)
        if "R4" in pmt[i]:
            data4.append(datan[i])
            data4.append(datan[i]+18)
            ws4.append(float(1./30.))
            ws4.append(float(1./30.))
            for n in range(1,5):
                data4.append(datan[i]+2*n)
                ws4.append(float(n+1)/30.)
                data4.append(datan[i]+(18-2*n))
                ws4.append(float(n+1)/30.)
        if "R5" in pmt[i]:
            data5.append(datan[i])
            data5.append(datan[i]+18)
            ws5.append(float(1./30.))
            ws5.append(float(1./30.))
            for n in range(1,5):
                data5.append(datan[i]+2*n)
                ws5.append(float(n+1)/30.)
                data5.append(datan[i]+(18-2*n))
                ws5.append(float(n+1)/30.)

    xhist1,xedges1 = histogram(data1,bins=nbins,weights=ws1,range=((xmin),(xmax)))
    xhist2,xedges2 = histogram(data2,bins=nbins,weights=ws2,range=((xmin),(xmax)))
    xhist3,xedges3 = histogram(data3,bins=nbins,weights=ws3,range=((xmin),(xmax)))
    xhist4,xedges4 = histogram(data4,bins=nbins,weights=ws4,range=((xmin),(xmax)))
    xhist5,xedges5 = histogram(data5,bins=nbins,weights=ws5,range=((xmin),(xmax)))

    x1 = array( [ (xedges1[i]+xedges1[i+1])/2. for i in range(len(xedges1)-1) ] )
    x2 = array( [ (xedges2[i]+xedges2[i+1])/2. for i in range(len(xedges2)-1) ] )
    x3 = array( [ (xedges3[i]+xedges3[i+1])/2. for i in range(len(xedges3)-1) ] )
    x4 = array( [ (xedges4[i]+xedges4[i+1])/2. for i in range(len(xedges4)-1) ] )
    x5 = array( [ (xedges5[i]+xedges5[i+1])/2. for i in range(len(xedges5)-1) ] )
    pl.errorbar(x1,xhist1,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x2,xhist2,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x3,xhist3,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x4,xhist4,yerr=0,drawstyle='steps-mid',capsize=0)
    pl.errorbar(x5,xhist5,yerr=0,drawstyle='steps-mid',capsize=0)
    
    net1 = sum(xhist1)
    net2 = sum(xhist2)
    net3 = sum(xhist3)
    net4 = sum(xhist4)
    net5 = sum(xhist5)
    rat5 = net1/net5
    rat4 = net2/net5
    rat3 = net3/net5
    rat2 = net4/net5
    pl.ylabel(yval,fontsize=20)
    pl.title(title,fontsize=24)

    textstr = '$N_R1$=%5.0f \n $N_R2$=%5.0f \n $N_R3$=%5.0f \n $N_R4$=%5.0f \n $N_R5$=%5.0f'% (net5,net4,net3,net2,net1)
    props=dict(boxstyle='round',facecolor='white',alpha=0.5)
    pl.text(0.75,0.95,textstr,fontsize=16,verticalalignment='top',transform=ax.transAxes,bbox=props)
    f.write('\n 5/1 = %3.3f \n 4/1 = %3.3f \n 3/1 = %3.3f \n 2/1 = %3.3f \n'% (rat5,rat4,rat3,rat2))
    f.close()


#function to create a waveform and use it to get events according to a definition (currently an old one; may need updating). Takes in double[] data=PE times and int[] coated=0 for uncoated, 1 for coated. 
#    returns double[] evtpes=an array of PE values for events, double[] stimes=array of starting times for those events with respect to the original event start time (for singlet/triplet), double[] evtpesun=the PE value for only uncoated PMTs for the events, and doublt[] coatratios= the uncoated/coated ratio for each event.
def getevents(data,coated):
    shift = np.random.uniform(-16000,16000)
    xmin=0
    xmax=16000
    nbins=np.linspace(xmin,xmax,8001)
    nevents = len(data)
    data2=[]
    data3=[]
    ws=[]
    ws3=[]
    data2 = data[:]

    ws = [1./30.]*len(data2)
    for i in range(0,nevents):
        data2[i] = data2[i]+shift
        data2.append(data2[i]+18)
        ws.append(float(1./30.))
        for n in range(1,5):
            nv = n+0.001-0.001
            data2.append(data2[i]+2*n)
            ws.append(float(n+1)/30.)
            data2.append(data2[i]+(18-2*n))
            ws.append(float(n+1)/30.)

        if (coated[i]==0):
            data3.append(data2[i])
            data3.append(data2[i]+18)
            ws3.append(float(1./30.))
            ws3.append(float(1./30.))
            for n in range(1,5):
                data3.append(data2[i]+2*n)
                ws3.append(float(n+1)/30.)
                data3.append(data2[i]+(18-2*n))
                ws3.append(float(n+1)/30.)
    #print(times,ws)                                            

    xhist,xedges = histogram(data2,bins=nbins,weights=ws,range=((xmin),(xmax)))
    xhist2,xedges2 = histogram(data3,bins=nbins,weights=ws3,range=((xmin),(xmax)))

    eventcount = 0
    evtpes = []
    evtpesun = []
    stimes = []
    evt = False
    ebin = -100
    for i in range(0,8000):
        if (xhist[i] > 0.4) and not evt:
            evt = True
            notstart = True
            notend = True
            sbin = i
            ebin = i
            holdpe = 0.
            holdun = 0.
            while (notstart):
                sbin = sbin -1
                if (sbin == 0):
                    sbin = 0
                    notstart = False
                elif (xhist[sbin] < .02):
                    notstart = False
            while (notend):
                ebin = ebin + 1
                if (ebin >=3996):
                    ebin = 4000
                    notend = False
                else:
                    next20 = xhist[ebin]+xhist[ebin+1]+xhist[ebin+2]+xhist[ebin+3]+xhist[ebin+4]
                    if (next20 < .02):
                        notend = False
            if (ebin-sbin)<45:
                lbin = ebin
            else:
                lbin = sbin + 45
            for n in range(sbin,lbin):
                holdpe = holdpe+xhist[n]
                holdun = holdun+xhist2[n]

            if holdpe < 3:
                evt = False
                ebin = 0
            else:
                holdtm = sbin*2-shift
                evtpes.append(holdpe)
                evtpesun.append(holdun)
                stimes.append(holdtm)
                if (i > ebin + 50) and evt and (xhist[i] < 0.02):
                    evt = False

    coatratios = []
    testval = 0.
    for n in range(0,len(evtpes)):
        testval = (evtpes[n]-evtpesun[n])/96.
        if testval <= 0.01:
                coatratios.append(100)
        else:
            testval = evtpesun[n]/testval
            testval = testval / 24.
            coatratios.append(testval)

    return evtpes,stimes,evtpesun,coatratios

#function to plot the event PEs histogram obtained from the previous function. Creates a figure on the given canvas (fig) in the given subplot (sub=1 or 2) subplot 1 creates linear y plots, subplot 2 log y. 
def plotPEs(evtpes,title,nbins,fig,sub,label):
    # now plot those evtpes                                                                                
    ax = fig.add_subplot(2,1,sub)

    pemin = min(evtpes)
    pemax = max(evtpes)
    thresh = pemax/3*2
    nbins = int(thresh/2)
    if nbins < 20:
        nbins = 50
    pehist,peedges = histogram(evtpes,bins=nbins,range=((pemin),(pemax)))
    pe = array( [ (peedges[i]+peedges[i+1])/2. for i in range(len(peedges)-1) ] )
    pl.errorbar(pe,pehist,yerr=0,marker='.',drawstyle='steps-mid',capsize=0)
    pl.xlabel('Count',fontsize=20)
    pl.title(title,fontsize=24)
    three7 = 0
    ten50 = 0
    inmax = 0
    for n in range(0,nbins):
        if peedges[n] > thresh:
            inmax = inmax + pehist[n]
        elif peedges[n] > 50:
            continue
        elif peedges[n] > 10:
            ten50 = ten50 + pehist[n]
        elif peedges[n] > 7:
            continue
        elif peedges[n] > 3:
            three7 = three7 + pehist[n]
        else:
            continue

    nevents = len(evtpes)
    mu=np.mean(evtpes)
    textstr = '$N$=%5.0f \n $\mu$=%.3f \n $N 3<x<7$ = %5.0f \n $N 10<x<50$ = %5.0f \n $N>%3.3f$ = %5.0f \n\
'% (nevents,mu,three7,ten50,thresh,inmax)
    props=dict(boxstyle='round',facecolor='white',alpha=0.5)
    pl.text(label,0.95,textstr,fontsize=16,verticalalignment='top',transform=ax.transAxes,bbox=props)
    if (sub==2):
        ax.set_yscale('log')

#Plots the coated/uncoated ratios. Also takes in a calculated average to compare to the one it calculates. 
def plotRat(ratios,ws,title,calculated,fig,sub,label):
    # now plot those ratios                         
    ax = fig.add_subplot(2,1,sub)
                                                       
    pemin = 0
    pemax = 10
    nbins = 200

    testval = 0
    testtot = 0
    for n in range(0,len(ratios)):
        testval = ratios[n]*ws[n]
        testtot = testtot+testval
    avg = testtot/sum(ws)

    pehist,peedges = histogram(ratios,bins=nbins,weights=ws,range=((pemin),(pemax)))
    pe = array( [ (peedges[i]+peedges[i+1])/2. for i in range(len(peedges)-1) ] )
    pl.errorbar(pe,pehist,yerr=0,marker='.',drawstyle='steps-mid',capsize=0)
    pl.xlabel('Count',fontsize=20)
    pl.title(title,fontsize=24)
    print min(ratios),max(ratios)

    nevents = len(ratios)
    mu=np.average(ratios,weights=ws)
    textstr = '$N$=%5.0f \n $\mu$=%.3f \n avg=%3.3f \n whole=%3.3f'% (nevents,mu,avg,calculated)
    props=dict(boxstyle='round',facecolor='white',alpha=0.5)
    pl.text(label,0.95,textstr,fontsize=16,verticalalignment='top',transform=ax.transAxes,bbox=props)
    if (sub==2):
        ax.set_yscale('log')

#calculates the ratio from the PEs per event and the uncoated PEs for those same events. returns the average ratio.
def getratio(evtpes,evtpesun):
    totalpes = sum(evtpes)
    totalpun = sum(evtpesun)
    ratio = (totalpun/24)/((totalpes-totalpun)/96)
    return ratio

#PMT response function. Takes in the energy (currently minimally considered) and the angle of incidence (cos) and returns True (photon converted to PE) or false (photon failed to convert).
def pmtresponse(eneg,cos):
    testval = np.random.uniform(0.0,1.0)
    if (cos < 0):
        cos = 0
    prob = 1.0*(cos**(0.5))
    if eneg > 4.0:
        return False
    elif eneg > 2.0 and testval < prob:
        return True
    else:
        return False

#function for reading in a file and processing it. returns allevts=event PEs,allunps=event PEs for uncoated only,alltimes=times for all PEs,allcu=PMT names for all PEs,allratios=uncaoted/coated ratios for all events
def processfile(infile,pdf2,var):
    print infile
    f0 = open(infile,"r")
    contents=f0.readlines()[0: ]
    
    unpes = []
    times = []
    petimes = []
    coatUn = []
    alltimes = []
    allcu = []
    allevts = []
    allunps = []
    allptimes = []
    allratios = []
    edpos = []
    npes  = []
    countpe = 0
    
    countun = 0
    countco = 0
    
    linen = 0
    nc = 0
    nevents = 0
    hc = 1239.841984
    
    for line in contents:
        line = line.strip()
        linen = linen+1
        #if (nevents > 500):                                                                                   
         #  continue                                                                                          
        if "Start" in line:
            if (nevents > 0):
                npes.append(float(countpe))
                
            #print(times[0],times[100],times[200])                                                         
                for count in range(0,5):
                    pesforevt,petimes,unpes,ratios = getevents(times,coatUn)
                    for i in range(0,len(pesforevt)):
                        allevts.append(pesforevt[i])
                        allunps.append(unpes[i])
                        allptimes.append(petimes[i])
                        allratios.append(ratios[i])
                    countun = countun+sum(unpes)
                    countco = countco+sum(pesforevt)-(sum(unpes))

                if (nevents%100 == 0):
                    print nevents,sum(allevts),sum(allunps)
                if (nc < 3) and (countpe > 1):
                    fig = pl.figure(figsize=(36,4))
                    ax=fig.add_subplot(1,1,1)
                    title="Drawn Waveform Event "+str(nevents)+"\n nPEs:"+str(countpe)
                    yvals = "Num PEs"
                    shift = 1000+100*var
                    buildwindow(times,coatUn,title,yvals,shift)
                    pdf2.savefig(fig)
                    pl.close()
                    nc = nc+1
            countpe = 0
            unpes=[]
            times =[]
            coatUn=[]
            petimes = []
            nevents = nevents + 1
            
        if "#" in line:
            continue
        
        elif "PMT" in line:
            if (line.count('\t') == 2):
                part,peen,time = line.split('\t',2)
                opeen = float(peen)
                countpe = countpe + 1
                times.append(float(time))
                alltimes.append(float(time))
                coatUn.append(part)
                allcu.append(part)
            elif (line.count('\t') == 3):
                part,peen,time,cos = line.split('\t',3)
                opeen = float(peen)
                cost = np.abs(float(cos))
                #cost=1
                if pmtresponse(opeen,cost):
                    countpe = countpe+1
                    times.append(float(time))
                    alltimes.append(float(time))
                    coatUn.append(part)
                    allcu.append(part)

    return allevts,allunps,alltimes,allcu,allratios

#Another function for calculating uncoated/coated ratios, using PE times and PMT names
def calcRatios(data,pmt):
    nevents = len(data)
    nett = 0
    for i in range(0,nevents):
        if "coated" not in pmt[i]:
            nett = nett+1

    ucT = float(nett)/24/(float(nevents-nett)/96)

    return ucT
    
#Same as above, but considers those PMTs off which are supposed to be off (for CCM120)
def wecalcRatios(data,pmt):
    nevents = len(data)
    nett = 0
    for i in range(0,nevents):
        if "coated" not in pmt[i]:
            nett = nett+1

    ucT = float(nett)/21/(float(nevents-nett)/78)

    return ucT
    
#Calculates the row ratios without producing a plot. Useful for saving time if you don't need visuals
def rowRatios(data,pmt,title,fn):
    filename = ["LaserRatios.txt","rerunLaserRatios.txt","fastRatios.txt","unweightFastRatios.txt"]
    f = open(filename[fn],"a")
    f.write(title)

    nevents = len(data)
    net1=0
    net2=0
    net3=0
    net4=0
    net5=0

    for i in range(0,nevents):
        if "R1" in pmt[i]:
            net1 = net1+1
        elif "R2" in pmt[i]:
            net2=net2+1
        elif "R3" in pmt[i]:
            net3=net3+1
        elif "R4" in pmt[i]:
            net4=net4+1
        elif "R5" in pmt[i]:
            net5=net5+1
       
    rat5 = float(net1)/net5
    rat4 = float(net2)/net5
    rat3 = float(net3)/net5
    rat2 = float(net4)/net5

    f.write('\n 5/1 = %3.3f \n 4/1 = %3.3f \n 3/1 = %3.3f \n 2/1 = %3.3f \n'% (rat5,rat4,rat3,rat2))
    f.close()

    return [rat5,rat4,rat3,rat2]

#process file function, but with some efficiency increases and fewer outputs. Not conducive for plot creation. Returns only PE times and PMT names. Also considers the SPE values and assigns an efficiency correction based on them (including off PMTs as off). Currently doesn't work (I forgot the file; will see about recreating it).
def fastprocessfile(infile):
    print infile
    alltimes = []
    allcu = []
    countpe = 0
    evtnum = 0

    #weights = speWeights("spevalues.txt")
    
    with open(infile,"r") as contents:
        for line in contents:
            line = line.strip()
            if "Start" in line:
                evtnum = evtnum+1
                if evtnum%1000 == 0:
                    print evtnum, countpe
#            if (evtnum > 5000):
#                continue
                    
            if "#" in line:
                continue
            
            elif "PMT" in line:
                if (line.count('\t') == 2):
                    part,peen,time = line.split('\t',2)
                    opeen = float(peen)
                    countpe = countpe + 1
                    alltimes.append(float(time))
                    allcu.append(part)
                elif (line.count('\t') == 3):
                    part,peen,time,cos = line.split('\t',3)
                    opeen = float(peen)
                    cost = np.abs(float(cos))
               #     we = returnweight(part,weights)
                    cost = cost*1
                #cost=1
                    if pmtresponse(opeen,cost):
                        countpe = countpe+1
                        alltimes.append(float(time))
                        allcu.append(part)
                
    return alltimes,allcu

#Same as above, but without the SPE weights (so exactly the same currently).
def unweightprocessfile(infile):
    print infile
    alltimes = []
    allcu = []
    countpe = 0
    evtnum = 0

    with open(infile,"r") as contents:
        for line in contents:
            line = line.strip()
            if "Start" in line:
                evtnum = evtnum+1
                if evtnum%1000 == 0:
                    print evtnum, countpe
                    
            if "#" in line:
                continue
            
            elif "PMT" in line:
                if (line.count('\t') == 2):
                    part,peen,time = line.split('\t',2)
                    opeen = float(peen)
                    countpe = countpe + 1
                    alltimes.append(float(time))
                    allcu.append(part)
                elif (line.count('\t') == 3):
                    part,peen,time,cos = line.split('\t',3)
                    opeen = float(peen)
                    cost = np.abs(float(cos))
                #cost=1
                    if pmtresponse(opeen,cost):
                        countpe = countpe+1
                        alltimes.append(float(time))
                        allcu.append(part)
                
    return alltimes,allcu

#Calculates the figures of Merit used for comparing the output of the Laser Simulation row and uncoated/coated ratios to the data values (listed below). 
def figureOfMerit(rowratios,ucratios,weighted):
    rowdata = [1.11, 1.06, 1.05, 0.954, 0.858, 0.909, 0.958, 0.904, 0.81, 0.884, 0.946, 0.902, 1.22, 1.16, 1.14, 1.03, 0.99, 1.06, 1.07, 1.00, 0.879, 1.02, 1.04, 1.00]
    ucdata = [0.98, 0.99, 0.98, 1.01, 1.06, 1.09]
    if (weighted) :
        rowdata = [1.11, 1.28, 1.05, 0.954, 0.858, 1.10, 0.958, 0.904, 0.81, 1.07, 0.946, 0.902, 1.22, 1.40, 1.14, 1.03, 0.99, 1.28, 1.07, 1.00, 0.879, 1.23, 1.04, 1.00]

    simplefig = 0
    squarefig = 0
    nonabsfig = 0

    if len(rowratios) != len(rowdata):
        print "Row ratios incomplete"
    elif len(ucdata) != len(ucratios):
        print "Coated/Uncoated ratios incomplete"
    else:
        for rr in range(0,len(rowdata)):
            simplefig = simplefig + abs(rowdata[rr]-rowratios[rr])/rowdata[rr]
            squarefig = squarefig + ((rowdata[rr]-rowratios[rr])**2)/rowdata[rr]
            nonabsfig = nonabsfig + (rowdata[rr]-rowratios[rr])/rowdata[rr]

        for uu in range(0,len(ucdata)):
            simplefig = simplefig + abs(ucdata[uu]-ucratios[uu])/ucdata[uu]
            squarefig = squarefig + ((ucdata[uu]-ucratios[uu])**2)/ucdata[uu]
            nonabsfig = nonabsfig + (ucdata[uu]-ucratios[uu])/ucdata[uu]

    return simplefig,squarefig,nonabsfig

#Weights the PMTs according to files that don't currently exist. Not very helpful until I correct that. 
def weightPMTs(nfiles):
    infile = ["213Top.txt", "213Middle.txt", "213Bottom.txt", "532Top.txt", "532Middle.txt", "532Bottom.txt"]
    holdvals = [""]*25
    avg = [0]*5
    for nn in range(0,5):
        avg[nn] = [0.]*24
    for n in range(0,nfiles):
        with open(infile[n],"r") as f0:
            linen = -1
            for line in f0:
                linen = linen+1
                holdvals = line.split(',',24)
                for i in range(0,24):
                    avg[linen][i]=avg[linen][i]+float(holdvals[i])

    mean = float(0)
    for n in range(0,5):
        for i in range(0,24):
            mean = mean+(avg[n][i])/120.
    for n in range(0,5):
        for i in range(0,24):
            avg[n][i] = avg[n][i]/mean/5

    return avg

#Same as above, except with the SPE values from the LED calibration. Also doesn't work without the source file.
def speWeights(infile):
    infile = "spevalues.txt"
    weis = [0]*5
    for nn in range(0,5):
        weis[nn] = [0.]*24
    with open(infile,"r") as f0:
        for line in f0:
            col,row,speval,pwr = line.split('\t',3)
            #print col,row,float(speval),pwr
            r = int(row)-1
            #c = int(col)-1
            c = int(col)
            if c==24:
                c=0
            if (int(pwr) == 0):
                weis[r][c] = 0
            else:
                spev = float(speval)
                errs = (spev)**(0.5)
                stdv = (spev-5)/(errs)
                weis[r][c] = norm.cdf(stdv)

    return weis
    
#calculates the PMT weight obtained from the SPE valuation above according to Row and column. Not currently useful.
def returnweight(pmtn,weights):
    r = pmtn.find('R')
    c = pmtn.find('C')
    row = 5-int(pmtn[r+1:r+2])
    column = int(pmtn[c+1:r])-1
    
    val = weights[row][column]
    return val

    
