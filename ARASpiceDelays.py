from __future__ import print_function
from ROOT import TCanvas, TGraph
from ROOT import gROOT
import ROOT
import os
import matplotlib.pyplot as plt
from scipy import signal,stats
from scipy.stats import norm
from scipy.fftpack import rfft, irfft, fftfreq
from scipy.signal import iirfilter,lfilter,butter,hilbert
from scipy import signal
from scipy.fftpack import fft
from scipy import optimize
from scipy.misc import derivative
from scipy.interpolate import interp1d, Akima1DInterpolator
import numpy as np
import sys
import math
from math import sin
from array import array
from pynverse import inversefunc
#from TimingCalibration import partial_derivative
#from VoltageCalibration import BlockCorrector, MeanEachBlock
#from AutomaticLoadData import LoadDataFromWeb, DownloadPedFile,LoadSpiceData
#sys.path.append('../')
from AutomaticLoadData import LoadSpiceData
import matplotlib
import scipy
import itertools
import warnings
import datetime
from scipy import interpolate
import ctypes
from scipy.interpolate import RegularGridInterpolator

warnings.simplefilter(action='ignore', category=FutureWarning)
ROOT.gSystem.Load("/home/kahughes/AraSim/libAra.so")
def R(cmd):
  return ROOT.gInterpreter.ProcessLine(cmd)
#ROOT.gInterpreter.AddIncludePath("/home/kahughes/AraSim")

def LoadModel(n0,nf,l):
    R('#include "/home/kahughes/AraSim/RayTrace.h"')
    R('#include "/home/kahughes/AraSim/RayTrace_IceModels.h"')
    R('#include "/home/kahughes/AraSim/Vector.h"')

    # attenuation model

    R('auto atten_model = boost::shared_ptr<basicAttenuationModel>( new basicAttenuationModel );');
    #exponential model
    R('auto refr_model = boost::shared_ptr<exponentialRefractiveIndex>(new exponentialRefractiveIndex('+str(n0)+','+str(nf)+','+str(l)+'));') #1.249,1.774,0.0163 #(1.353,1.78,0.0160)


    R('RayTrace::TraceFinder tf(refr_model, atten_model);')
    R('Vector src; Vector rec; std::vector<RayTrace::TraceRecord> paths;')
def trace(r_src, z_src, z_rec=-200):
    ROOT.src.SetXYZ(r_src,0,z_src);
    ROOT.rec.SetXYZ(0,0,z_rec);
    sol_cnt = ctypes.c_int()
    sol_err = ctypes.c_int()
    ROOT.paths = ROOT.tf.findPaths(ROOT.src, ROOT.rec, 0.5, ROOT.TMath.Pi()/2, sol_cnt, sol_err, ROOT.RayTrace.NoReflection,0.001)

    times = []
    for path in ROOT.paths:
        times.append({'tof':path.pathTime*1e9,'dist': path.pathLen, 'atten': path.attenuation})
    proptimes = times[0]["tof"] if len(times) else -1
    return proptimes

def LoadSpice():
    f = ROOT.TFile("Sp18Log.root")
    myTree = f.Get("N")
    #plt.figure(1,facecolor='w')
    times = []
    depths = []
    for entry in myTree:
        utcday = entry.utcday
        utchr = entry.utchr
        utcmin = entry.min
        #print(utcday,utchr,utcmin)
        date = datetime.datetime(2018,12,int(utcday)-334,int(utchr),int(utcmin),0)
        unix = (date-datetime.datetime(1970,1,1)).total_seconds()
        times.append(unix)
        #print(date)
        #print(unix)
        depth = entry.z
        depths.append(depth)
        #times.append(utcmin)
        #depths.append(depth)
        #print(utcday,utchr,utcmin)
        #plt.scatter(utcmin,depth)
    #plt.show()

    myTree2 = f.Get("T")
    for entry in myTree2:
        unixtime = entry.unixtime-46800+86400
        #print(unixtime)
        #times.append(unixtime)
        utcz = entry.utcz
        #depths.append(utcz)
    return(times,depths)

def my_chi2(s_depths,s_times,my_depths,my_delays):
    interp_times = interpolate.interp1d(s_depths,s_times)
    these_delays = my_delays[my_depths>-1400]
    these_depths = my_depths[my_depths>-1400]
    these_delays = these_delays[these_depths<-1000]
    these_depths = these_depths[these_depths<-1000]
    #print(interp_times(-1100))
    diff = np.abs(s_times[0]-my_delays[0])
    cable_delays = np.linspace(diff-5,diff+5,10)
    chi_all = []
    #print(my_depths)

    these_delays = these_delays[0::10]
    these_depths = these_depths[0::10]
    #for d in cable_delays:
    chi=0
    for d in these_depths:
        #print(d)
        expected_t = interp_times(d)
        #print(np.where(these_depths==d)[0])
        observed_t = these_delays[np.where(these_depths==d)[0]]
        #print(expected_t,observed_t)
        chi=chi+((observed_t-expected_t)**2)
    chi_all=chi
    """
    for c in cable_delays:
        chi = 0
        #print(len(my_depths))
        #plt.figure(1)
        #plt.scatter(my_delays,my_depths)
        #plt.plot(s_times,s_depths)
        #plt.show()
        for d in these_depths:
            #print(d)
            expected_t = interp_times(d)-c
            #print(np.where(these_depths==d)[0])
            observed_t = these_delays[np.where(these_depths==d)[0]]
            #print(expected_t,observed_t)
            chi=chi+((observed_t-expected_t)**2)
        #print(chi)
        chi_all.append(chi[0])
    """
    #plt.figure(0)
    #plt.plot(cable_delays,chi_all)
    #plt.show()
    return(np.min(chi_all))
def smooth(y,box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y,box,mode='same')
    return(y_smooth)

def AverageSpice(delays, depths):
    this_mean = np.mean(delays)
    my_vals = np.where(np.abs(delays-this_mean)<5)
    depths = depths[my_vals]
    delays = delays[my_vals]

    jumpsize = 25

    bin_delay = np.zeros(len(delays)/jumpsize)
    bin_depth = np.zeros(len(delays)/jumpsize)
    #print(len(delay)/50)
    for l in range(0,len(delays)/jumpsize):
        bin_delay[l]=np.mean(delays[l*jumpsize:l*jumpsize+jumpsize])
        bin_depth[l]=np.mean(depths[l*jumpsize:l*jumpsize+jumpsize])

    plt.figure(1,facecolor='w')
    plt.scatter(delays,depths,color='dodgerblue')
    #plt.plot(bin_delay,bin_depth,lw=2.5)
    plt.plot(smooth(delays,250),depths,lw=2.5)

    plt.show()


def SpiceChi(this_delay,this_depth):
    #PA_depth = {0:-172.635,1:-173.65,2:-174.66,3:-175.68,4:-176.70,5:-189.5,6:-178.75,7:-180.79}
    this_delay = np.asarray(this_delay)
    this_depth = np.asarray(this_depth)

    this_mean = np.mean(this_delay)
    my_vals = np.where(np.abs(this_delay-this_mean)<15)
    this_depth = this_depth[my_vals]
    this_delay = this_delay[my_vals]

    AV_delay = smooth(this_delay,150)
    #AverageSpice(this_delay,this_depth)
    plt.scatter(this_delay,this_depth,color='dodgerblue')
    plt.plot(AV_delay,this_depth)
    plt.show()
    LoadModel(1.326,1.78,0.0202)
    #delays = np.load('good_delays.npy')
    #depths = np.load('good_depths.npy')
    #this_delay = delays[4,:]
    #this_depth = depths[4,:]
    #print(this_delay,this_depth)
    distances = np.linspace(20.0,25.0,18)
    #A_depths = np.linspace(3.0,5.0,20)
    A_depths = np.linspace(-190.0,-192.0,20)
    #ARA_depths = np.asarray([-191.28])
    #SPICE_depths = np.linspace(-900,-1400,30)
    AV_delay = AV_delay[(this_depth<-1020)&(this_depth>-1450)]
    SPICE_depths = this_depth[(this_depth<-1020)&(this_depth>-1450)]
    SPICE_depths = SPICE_depths[0::100]
    AV_delay = AV_delay[0::100]
    SPICE_times = np.zeros(len(SPICE_depths))
    highval = [-3,-4,-8,-10,-170,-15,-25]
    lowval = [-6,-10,-15,-20,-195,-26,-35]
    #all_chis = np.zeros([len(A_depths),len(B_depths),len(distances)])
    all_chis = np.zeros([len(A_depths),len(distances)])
    #this_depth = this_depth[np.argwhere((this_delay<50) & (this_delay>20))]
    #this_delay = this_delay[np.argwhere((this_delay<50) & (this_delay>20))]
    best_chi = 100000.0
    print(len(SPICE_depths))

    r =  23.24
    #z1 = -195.00
    #z2 = -190.0
    counter = 0

    for z_s in SPICE_depths:
        #print('')
        t1 = trace(4140.0,z_s,-191+29.85)#ch 8 is the closer one
        t2 = trace(4140.0,z_s,-191)#ch 24 is farther
        SPICE_times[counter]=t2-t1

        if(z_s==SPICE_depths[10]):
            print('here')
            lineup=SPICE_times[10]-AV_delay[10]
            #lineup= -113.68
        #print(SPICE_times[counter],this_delay[np.where(z_s==SPICE_depths)[0][0]]+lineup)
        #all_chis[np.where(A_depths==z1)[0],np.where(B_depths==z2)[0],np.where(distances==r)[0]]=(all_chis[np.where(A_depths==z1)[0],np.where(B_depths==z2)[0],np.where(distances==r)[0]]+
        #(AV_delay[np.where(z_s==SPICE_depths)[0][0]]-SPICE_times[counter]+lineup)**2/(0.1**2+0.1**2))
        counter = counter +1

    plt.figure(1,facecolor='w')
    plt.scatter(this_delay+lineup,this_depth,color='dodgerblue')
    #plt.plot(AV_delay,SPICE_depths,lw=2.5,label='Rolling Average',color='darkblue')
    plt.plot(SPICE_times,SPICE_depths,lw=2.5,label='Best Fit',color='red')
    plt.legend()
    plt.xlabel('Time Delay (ns)')
    plt.ylabel('SPIceCore Depth (m)')
    plt.title('Cable Delay='+str(lineup))
    plt.show()
    cal_delays = {(8, 25): 310.7580988219895, (16, 17): 73.63767997382199, (1, 24): 46.942653795811516, (24, 16): 0.3178174083769634, (25, 16): -128.64406086387436, (0, 25): 315.0504744764398, (1, 17): 120.86690117801047, (8, 9): 136.76464332460733, (9, 16): 45.349394633507856, (0, 16): 186.42203861256544, (1, 9): 2.004826570680628, (8, 17): 255.68921793193718, (0, 8): 4.308000654450262, (9, 25): 173.9465804973822, (25, 17): -55.084505890052355, (0, 1): 139.06781740837695, (24, 25): 128.9462532722513, (8, 24): 181.76497054973822, (1, 25): 175.93578206806282, (1, 16): 47.369846204188484, (0, 24): 186.0885962041885, (0, 17): 259.9815935863874, (24, 17): 73.87737238219896, (1, 8): -134.8066917539267, (9, 17): 118.90894960732984, (0, 9): 141.0882689790576, (9, 24): 44.98470222513089, (8, 16): 182.09841295811518}
    for r in distances:
        print('next r is', r)
        for z1 in A_depths:

            #for z2 in B_depths:
            #print(z2)

            counter = 0
            #print(r)
            for z_s in SPICE_depths:
                #print('')
                my_zval = -194.5
                t1 = trace(4160.0,z_s,my_zval)#ch 8 is the closer one
                t2 = trace(4160.0,z_s,my_zval+29.85)#ch 24 is farther
                SPICE_times[counter]=t2-t1
                if(z_s==SPICE_depths[0]):
                    #print('here')
                    lineup=SPICE_times[0]-AV_delay[0]
                #print(SPICE_times[counter],this_delay[np.where(z_s==SPICE_depths)[0][0]]+lineup)
                all_chis[np.where(A_depths==z1)[0],np.where(distances==r)[0]]=(all_chis[np.where(A_depths==z1)[0],np.where(distances==r)[0]]+
                (AV_delay[np.where(z_s==SPICE_depths)[0][0]]-SPICE_times[counter]+lineup)**2/(0.1**2+0.1**2))
                counter = counter +1
            #print(all_chis[np.where(A_depths==z1)[0],np.where(B_depths==z2)[0],np.where(distances==r)[0]])

            if(all_chis[np.where(A_depths==z1)[0],np.where(distances==r)[0]]<len(SPICE_depths)):
                print(all_chis[np.where(A_depths==z1)[0],np.where(distances==r)[0]])
                print(r,z1,z2)
                print(z2-z1)
                print('')
            if(all_chis[np.where(A_depths==z1)[0],np.where(distances==r)[0]]<best_chi):

                best_chi = all_chis[np.where(A_depths==z1)[0],np.where(distances==r)[0]]

            #if(r>4155 and z<-191):

            #plt.figure(0,facecolor='w')
            #plt.scatter(this_delay[this_depth>-1400],this_depth[this_depth>-1400],color='dodgerblue')
            #plt.plot(SPICE_times[SPICE_depths<-1000]-SPICE_times[15],SPICE_depths[SPICE_depths<-1000],lw=2.0)
            #plt.xlabel('Time delay (ns)')
            #plt.ylabel('SPIceCore Depth (m)')
            #plt.title('Original')
            #plt.show()
            #6.8,0.9,4.0,3.2

            #print(this_delay)
            #print(this_depth)
            #all_chis[np.where(ARA_depths==z)[0],np.where(distances==r)[0]]=my_chi2(SPICE_depths,SPICE_times,this_depth,this_delay)
            #plt.figure(0,facecolor='w')
            #plt.plot(SPICE_times,SPICE_depths)
            #plt.scatter(this_delay,this_depth)
            #plt.show()
        #plt.show()
    min_chi2=np.unravel_index(all_chis.argmin(), all_chis.shape)
    print('best values are', distances[min_chi2[1]],A_depths[min_chi2[0]],np.min(all_chis))
    plt.figure(0)
    plt.imshow(all_chis,aspect='auto',origin='upper',interpolation='None',extent=[np.min(distances),np.max(distances),np.min(A_depths),np.max(A_depths)])#,vmin=0,vmax=21)
    #plt.plot(distances,all_chis[0])
    plt.show()
    #plt.figure(1)
    #plt.imshow(all_chis[:,0,:],aspect='auto',origin='upper',interpolation='None',extent=[np.min(distances),np.max(distances),np.min(A_depths),np.max(A_depths)])#,vmin=0,vmax=21)


    #plt.show()

def FindCalDelays():
    times,depths = LoadSpice()
    f_depths = interpolate.interp1d(times,depths)
    event_times = np.load("cAL_event_times.npy")
    event_depths = []
    event_delays = []
    max_cors = []
    #cutoffs = {0:260.0,1:175.0,8:310.0,9:210,16:475.0,17:375.0,24:450.0,25:350.0}
    for i in range(2500,len(t0[:,0])):#d.N()):
        print(i)


        if(volt0[i,-1]==0.0 or volt1[i,-1]==0.0):
            continue
        t_cal_even = np.arange(t0[i,0],t0[i,-1]-1.0,0.3125)
        f = Akima1DInterpolator(t0[i,:],volt0[i,:])
        this_v0= f(t_cal_even)

        f = Akima1DInterpolator(t1[i,:],volt1[i,:])
        this_v1 = f(t_cal_even)

        this_v0 = signal.resample(this_v0,len(this_v0)*10)
        this_v1,t_cal_even = signal.resample(this_v1,len(this_v1)*10,t=t_cal_even)
        #if(np.max(this_v0)<150):
        #    continue
        dt = t_cal_even[1]-t_cal_even[0]

        this_v0[int(cutoffs[ch0]/dt):] = 0.0
        this_v1[int(cutoffs[ch1]/dt):] = 0.0

        cor = signal.correlate(this_v0/np.std(this_v0),this_v1/np.std(this_v1))/len(this_v0)
        max_cors.append(np.max(cor))
        #if(np.max(cor)>0.65):
        this_delay = (np.argmax(cor)-np.size(cor)/2.)*dt
        event_delays.append(this_delay)
        event_depths.append(f_depths(event_times[i]))

    return(event_delays,event_depths)

def FindSPICEDelays(t0,t1,volt0,volt1,ch0,ch1):
    times,depths = LoadSpice()
    
    f_depths = interpolate.interp1d(times,depths)
    event_times = np.load("/home/kahughes/PA_Analysis/data/SPICE_event_times.npy")
    event_depths = []
    event_delays = []
    max_cors = []
    cutoffs = {0:280.0,1:195.0,8:330.0,9:230,16:495.0,17:395.0,24:470.0,25:370.0}
    for i in range(2500,len(t0[:,0])):#d.N()):
        #print(i)


        if(volt0[i,-1]==0.0 or volt1[i,-1]==0.0):
            continue
        t_cal_even = np.arange(t0[i,0]+1.0,t0[i,-1]-1.0,0.03125)
        f = Akima1DInterpolator(t0[i,:],volt0[i,:])
        this_v0= f(t_cal_even)
        #t1[i,:]=t1[i,:]-t0[i,0]
        f = Akima1DInterpolator(t1[i,:],volt1[i,:])
        this_v1 = f(t_cal_even)

        #this_v0 = signal.resample(this_v0,len(this_v0)*10)
        #this_v1,t_cal_even = signal.resample(this_v1,len(this_v1)*10,t=t_cal_even)
        #if(np.max(this_v0)<150):
        #    continue
        dt = t_cal_even[1]-t_cal_even[0]

        this_v0[int(cutoffs[ch0]/dt):] = 0.0
        this_v1[int(cutoffs[ch1]/dt):] = 0.0
        """
        if(ch0==25 or ch1==25):
            #print(this_v1)
            #print(t0[i,:])
            #print(t1[i,:])
            #print(t_cal_even)
            plt.plot(t_cal_even,this_v0)
            plt.plot(t_cal_even,this_v1)
            plt.show()
            plt.plot(t0[i,:],volt0[i,:])
            plt.plot(t1[i,:],volt1[i,:])
            plt.show()
        """
        #plt.figure(1)
        #plt.plot(t_cal_even,this_v0)
        #plt.plot(t_cal_even,this_v1)
        #plt.show()

        cor = signal.correlate(this_v0/np.std(this_v0),this_v1/np.std(this_v1))/len(this_v0)
        max_cors.append(np.max(cor))
        #if(ch1==16):
        #    plt.figure(1)
        #    plt.plot(t_cal_even,this_v0)
        #    plt.plot(t_cal_even,this_v1)
        #    plt.show()
        #    plt.plot(np.linspace(0,len(cor),len(cor)),cor)
        #    plt.show()
        #if(np.max(cor)>0.42 or i>len(t0[:,0])*.80):
        this_delay = (np.argmax(cor)-np.size(cor)/2.)*dt
        event_delays.append(this_delay)
        event_depths.append(f_depths(event_times[i]))

    #plt.figure(0,facecolor='w')
    #plt.hist(max_cors,bins=30)
    #colors={0:'dodgerblue',8:'springgreen',16:'red',24:'lavender'}
    #plt.scatter(np.asarray(event_delays)-113.68,event_depths,color=colors[ch0],label='Ch'+str(ch0)+' and Ch'+str(ch1))
    #plt.ylim([-1500,0.0])
    #plt.xlabel('Time Delays')
    #plt.ylabel('SPIceCore Depth')

    #plt.figure(1)
    #plt.hist(max_cors,bins=30)

    #plt.show()



    #all_cor_maxes = np.asarray(all_cor_maxes)
    #np.save('ARA_good_delays.npy',delays)
    #np.save('ARA_good_depths.npy',depths)
    return(event_delays,event_depths)
def TwoStringSpiceChi(delays,depths,pairs,ch1,ch2):


    SPICE_depths = np.linspace(-1050,-1450,20)
    interp_delays = {}
    for i in pairs:
        this_delay = np.asarray(delays[i])
        this_depth = np.asarray(depths[i])
        this_mean = np.mean(this_delay)
        my_vals = np.where(np.abs(this_delay-this_mean)<15)
        this_depth = this_depth[my_vals]
        this_delay = this_delay[my_vals]

        AV_delay = smooth(this_delay,350)
        f_int = interpolate.interp1d(this_depth,AV_delay)
        temp_delay = np.zeros(len(SPICE_depths))
        for d in range(0,len(SPICE_depths)):
            temp_delay[d]=f_int(SPICE_depths[d])
        interp_delays.update({i:np.asarray(temp_delay)})

        #plt.figure(0,facecolor='w')
        #plt.scatter(this_delay,this_depth,color='dodgerblue')
        #plt.plot(temp_delay,SPICE_depths,color='darkblue',lw=2.0)
        #plt.show()


    z_guesses = {0:-195.31,8:-194.52,24:-191.28,16:-177.72}
    cable_delays = {8:-114.78,24:-117.08}

    distances = np.linspace(22.0,24.0,10)
    #depths = np.linspace(-190.0,-195.0,10)
    #depths1 = np.linspace(z_guesses[ch1]-0.5,z_guesses[ch1]+1.5,10)
    #depths2 = np.linspace(z_guesses[ch2]-0.5,z_guesses[ch2]+1.5,10)
    #depth_diff = np.linspace(4.0,5.5,30)
    depth_diff = np.asarray([4.7])
    cables = np.linspace(0.0,7.0,10)
    #cables = np.asarray([-0.43])
    print(pairs)
    #model_times_b = np.zeros([len(distances),len(depths),len(SPICE_depths)])
    #model_times_t1 = np.zeros([len(distances),len(depths),len(SPICE_depths)])
    #model_times_t2 = np.zeros([len(distances),len(depths),len(SPICE_depths)])
    LoadModel(1.326,1.78,0.0202)

    P = -114.78
    Q = -117.08
    good_times = []
    counter = 0
    best_r = 20.0#21.33
    best_zd = 3.68#4.73
    best_c = 18.68#7.52
    for zs in SPICE_depths:
        t1 = trace(4160.0,zs,-164.67)#ch 8 is the closer one
        t2 = trace(4160+best_r,zs,-164.67+best_zd)#ch 24 is farther
        good_times.append(t2-t1+best_c)
        print(t2-t1+best_c,interp_delays[(9,25)][counter])
        counter = counter +1
    plt.figure(1,facecolor='w')
    this_delay = np.asarray(delays[9,25])
    this_depth = np.asarray(depths[(9,25)])
    plt.scatter(-this_delay,this_depth,color='dodgerblue')
    plt.plot(good_times,SPICE_depths,lw=2.0)
    plt.plot(-interp_delays[(9,25)],SPICE_depths,lw=2.0)

    plt.title(i)
    plt.show()

    rv = 24.05
    z1 = 6.0
    R = 0.0
    all_cables={(8,9):P,(24,25):Q,(8,24):R,(9,25):-P+R+Q,(8,25):R+Q,(9,24):R-P}
    solid_r = 4160.0
    solid_z = -164.67
    cal_delays = {(8, 25): 310.7580988219895, (16, 17): 73.63767997382199, (1, 24): 46.942653795811516, (24, 16): 0.3178174083769634, (25, 16): -128.64406086387436, (0, 25): 315.0504744764398, (1, 17): 120.86690117801047, (8, 9): 136.76464332460733, (9, 16): 45.349394633507856, (0, 16): 186.42203861256544, (1, 9): 2.004826570680628, (8, 17): 255.68921793193718, (0, 8): 4.308000654450262, (9, 25): 173.9465804973822, (25, 17): -55.084505890052355, (0, 1): 139.06781740837695, (24, 25): 128.9462532722513, (8, 24): 181.76497054973822, (1, 25): 175.93578206806282, (1, 16): 47.369846204188484, (0, 24): 186.0885962041885, (0, 17): 259.9815935863874, (24, 17): 73.87737238219896, (1, 8): -134.8066917539267, (9, 17): 118.90894960732984, (0, 9): 141.0882689790576, (9, 24): 44.98470222513089, (8, 16): 182.09841295811518}
    deep_delays ={(8, 25): -60.765625, (16, 17): 85.203125, (1, 24): -265.609375, (24, 16): -39.078125, (25, 16): -120.796875, (0, 25): -97.296875, (1, 17): -219.453125, (8, 9): 80.734375, (9, 16): -262.328125, (0, 16): -218.109375, (1, 9): -42.390625, (8, 17): -96.359375, (0, 8): -36.546875, (9, 25): -141.484375, (25, 17): -35.609375, (0, 1): 86.546875, (24, 25): 81.703125, (8, 24): -142.515625, (1, 25): -183.890625, (1, 16): -304.671875, (0, 24): -179.046875, (0, 17): -132.890625, (24, 17): 46.140625, (1, 8): -123.109375, (9, 17): -177.109375, (0, 9): 44.203125, (9, 24): -223.234375, (8, 16): -181.578125}

    best_chi = 100000
    dist = np.linspace(20.0,23.0,10)#max distance should be ~29.3
    zdep = np.linspace(0.0,5.0,20)
    cal_dist = np.linspace(30.0,32.0,10)
    mini_chi2 = np.zeros([len(zdep),len(dist),20])
    cr = 0

    eastings = [32355.0,32337.6,32420.0,32361.5,32290.0,32272.7,42358.93,44884.42]
    northings = [39734.5,39795.1,39745.5,39654.5,39710.6,39582.8,48974.22,51444.88]
    #plt.figure(0,facecolor='w')
    #plt.scatter(eastings,northings)
    #plt.show()
    """
    for r in dist:
        print(r,cr)
        cz = 0
        for z in zdep:
        #for cal_r in cal_dist:
            cc = 0
            for c in np.linspace(0.0,13.0,20):
                counter = 0
                for zs in SPICE_depths:

                    t1 = trace(solid_r,zs,solid_z)#ch 8 is the closer one
                    t2 = trace(solid_r+r,zs,solid_z+z)#ch 24 is farther

                    mini_chi2[cz,cr,cc]=mini_chi2[cz,cr,cc]+(t2-t1+c+interp_delays[(9,25)][counter])**2
                    #counter = 0
                    counter=counter +1

                #t1 = trace(38.5,-173.28,solid_z+z)
                #t2 = trace(38.5+29.0,-173.28,solid_z)
                #mini_chi2[cz,cr,cc]=mini_chi2[cz,cr,cc]+(t2-t1+c-cal_delays[(9,25)])**2

                #t1 = trace(5240.0,-1370.0,solid_z)
                #t2 = trace(5240.0+r,-1370.0,solid_z+z)
                #mini_chi2[cz,cr,cc]=mini_chi2[cz,cr,cc]+(t2-t1+c+deep_delays[(9,25)])**2
                #print(t2-t1+c,deep_delays[(9,25)])

                #print('')
                if(mini_chi2[cz,cr,cc]<1.0):
                    print(r,z,c,mini_chi2[cz,cr,cc])
                    best_chi=mini_chi2[cz,cr,cc]
                cc = cc+1

            cz = cz+1
        cr = cr+1
    """
    """
    plt.figure(0,facecolor='w')
    plt.imshow(mini_chi2,aspect='auto',origin='lower',interpolation='None',extent=[np.min(dist),np.max(dist),np.min(cal_dist),np.max(cal_dist)],vmin=0,vmax=120)

    plt.colorbar()
    plt.xlabel('Difference in r (m) between strings')
    plt.ylabel('Difference in z (m) between strings')
    #plt.axvline(distances[min_chi2[1]],color='white')
    #plt.axhline(depth_diff[min_chi2[0]],color='white')
    plt.show()
    """
    cal_delays = {(8, 25): 310.7580988219895, (16, 17): 73.63767997382199, (1, 24): 46.942653795811516, (24, 16): 0.3178174083769634, (25, 16): -128.64406086387436, (0, 25): 315.0504744764398, (1, 17): 120.86690117801047, (8, 9): 136.76464332460733, (9, 16): 45.349394633507856, (0, 16): 186.42203861256544, (1, 9): 2.004826570680628, (8, 17): 255.68921793193718, (0, 8): 4.308000654450262, (9, 25): 173.9465804973822, (25, 17): -55.084505890052355, (0, 1): 139.06781740837695, (24, 25): 128.9462532722513, (8, 24): 181.76497054973822, (1, 25): 175.93578206806282, (1, 16): 47.369846204188484, (0, 24): 186.0885962041885, (0, 17): 259.9815935863874, (24, 17): 73.87737238219896, (1, 8): -134.8066917539267, (9, 17): 118.90894960732984, (0, 9): 141.0882689790576, (9, 24): 44.98470222513089, (8, 16): 182.09841295811518}

    chi2 =np.zeros([len(cal_dist),len(distances),len(cables)])
    P = -114.78
    Q = -117.08
    best_chi = 1000000.0
    cr = 0
    z1 = 3.33
    for rv in distances:
        print(rv)
        cz1 = 0
        for cal_r in cal_dist:
            cc = 0
            for R in cables:
                #print(cr,cz1,cc)
                all_cables={(8,9):P,(24,25):Q,(8,24):R,(9,25):-P+R+Q,(8,25):R+Q,(9,24):R-P}
                #print(all_cables[(8,24)])

                for i in pairs:
                    #if(i==(9,25)):
                    #    continue
                    test_z = -194.6
                    test_r = 4160.0
                    zt1 = 29.85
                    zt2 = 29.94
                    if(i[0]==ch1):
                        this_z1 = test_z
                        this_r1 = test_r
                    if(i[0]==ch1+1):
                        this_z1 = test_z+zt1
                        this_r1 = test_r
                    if(i[0]==ch2):
                        this_z1 = test_z+z1
                        this_r1 = test_r+rv

                    if(i[1]==ch1+1):
                        this_z2 = test_z+zt1
                        this_r2 = test_r
                    if(i[1]==ch2):
                        this_z2 = test_z+z1
                        this_r2 = test_r+rv
                    if(i[1]==ch2+1):
                        this_z2 = test_z+z1+zt2
                        this_r2 = test_r+rv
                    counter = 0
                    #print(i,this_z1,this_r1,this_z2,this_r2)
                    for zs in SPICE_depths:

                        t1 = trace(this_r1,zs,this_z1)#ch 8 is the closer one
                        t2 = trace(this_r2,zs,this_z2)#ch 24 is farther
                        chi2[cz1,cr,cc]=chi2[cz1,cr,cc]+(t1-t2-all_cables[i]-interp_delays[i][counter])**2

                        #print(i,t1-t2-all_cables[i]-interp_delays[i][counter])
                        counter=counter +1

                    this_r1 = {24:38.5,25:38.5,8:38.5+cal_r,9:38.5+cal_r}
                    this_r2 = {24:38.5,25:38.5,8:38.5+cal_r,9:38.5+cal_r}
                    t1 = trace(this_r1[i[0]],-173.28,this_z1)
                    t2 = trace(this_r2[i[1]],-173.28,this_z2)
                    #print(t2,t1,all_cables[i],cal_delays[i])
                    #print(t1-t2-all_cables[i],cal_delays[i])
                    chi2[cz1,cr,cc]=chi2[cz1,cr,cc]+(t1-t2-all_cables[i]-cal_delays[i])**2

                #print(chi2[cz1,cr,cc])
                if(chi2[cz1,cr,cc]<best_chi):
                    print('new best!')
                    print(rv,cal_r,R,chi2[cz1,cr,cc])
                    best_chi=chi2[cz1,cr,cc]
                    """
                    for i in pairs:
                        test_z = -194.6
                        test_r = 4160.0
                        zt1 = 29.85
                        zt2 = 29.94
                        if(i[0]==ch1):
                            this_z1 = test_z
                            this_r1 = test_r
                        if(i[0]==ch1+1):
                            this_z1 = test_z+zt1
                            this_r1 = test_r
                        if(i[0]==ch2):
                            this_z1 = test_z+z1
                            this_r1 = test_r+rv

                        if(i[1]==ch1+1):
                            this_z2 = test_z+zt1
                            this_r2 = test_r
                        if(i[1]==ch2):
                            this_z2 = test_z+z1
                            this_r2 = test_r+rv
                        if(i[1]==ch2+1):
                            this_z2 = test_z+z1+zt2
                            this_r2 = test_r+rv
                        counter = 0
                        good_times = []
                        for zs in SPICE_depths:

                            t1 = trace(this_r1,zs,this_z1)#ch 8 is the closer one
                            t2 = trace(this_r2,zs,this_z2)#ch 24 is farther
                            good_times.append(t1-t2-all_cables[i])

                            chi2[cz1,cr,cc]=chi2[cz1,cr,cc]+(t1-t2-all_cables[i]-interp_delays[i][counter])**2/(0.3**2+0.1**2)

                            #print(i,t1-t2-all_cables[i]-interp_delays[i][counter])
                            counter=counter +1
                        plt.figure(1,facecolor='w')
                        plt.plot(good_times,SPICE_depths)
                        plt.plot(interp_delays[i],SPICE_depths)
                        plt.title(i)
                        plt.show()
                    """
                cc=cc+1

            cz1 = cz1+1
        cr=cr+1
    np.save('chi2_full.npy',chi2)
    min_chi2=np.unravel_index(chi2.argmin(), chi2.shape)
    print(min_chi2)
    #print('best values are', rvals[min_chi2[1]],zvals[min_chi2[0]],np.min(chi2))

    plt.figure(0,facecolor='w')
    #plt.imshow(chi2[:,:,0],aspect='auto',origin='lower',interpolation='None',extent=[np.min(distances),np.max(distances),np.min(depth_diff),np.max(depth_diff)],vmin=0,vmax=120)
    plt.imshow(chi2[0,:,:],aspect='auto',origin='lower',interpolation='None',extent=[np.min(cables),np.max(cables),np.min(distances),np.max(distances)],vmin=0,vmax=120)

    #plt.plot(distances,all_chis[0])
    plt.colorbar()
    plt.xlabel('Difference in r (m) between strings')
    plt.ylabel('Difference in z (m) between strings')
    plt.axvline(distances[min_chi2[1]],color='white')
    plt.axhline(depth_diff[min_chi2[0]],color='white')
    plt.show()

def TwoStringFitter(delays,depths,pairs,ch1,ch2):
    print('beginning two string fitter')
    cal_delays = {(8, 25): 310.7580988219895, (16, 17): 73.63767997382199, (1, 24): 46.942653795811516, (24, 16): 0.3178174083769634, (25, 16): -128.64406086387436, (0, 25): 315.0504744764398, (1, 17): 120.86690117801047, (8, 9): 136.76464332460733, (9, 16): 45.349394633507856, (0, 16): 186.42203861256544, (1, 9): 2.004826570680628, (8, 17): 255.68921793193718, (0, 8): 4.308000654450262, (9, 25): 173.9465804973822, (25, 17): -55.084505890052355, (0, 1): 139.06781740837695, (24, 25): 128.9462532722513, (8, 24): 181.76497054973822, (1, 25): 175.93578206806282, (1, 16): 47.369846204188484, (0, 24): 186.0885962041885, (0, 17): 259.9815935863874, (24, 17): 73.87737238219896, (1, 8): -134.8066917539267, (9, 17): 118.90894960732984, (0, 9): 141.0882689790576, (9, 24): 44.98470222513089, (8, 16): 182.09841295811518}
    """
    (8, 24): 181.76497054973822
    (9, 25): 173.9465804973822

    (8, 9): 136.76464332460733
    (24,25):128.9462532722513

    (9, 24): 44.98470222513089
    (8, 25): 310.7580988219895
    """
    interp = LoadPropTimes()
    SPICE_depths = np.linspace(-1050,-1450,20)
    interp_delays = {}
    for i in pairs:
        cut_vals = {(9,25):2.0,(8,25):15.0,(8,9):15.0,(8,24):2.0,(9,24):15.0,(24,25):15.0}
        linear_fit = {(9,25):1,(8,25):0,(8,9):0,(8,24):1,(9,24):0,(24,25):0}
        print(i)
        #i=(9,24)
        this_delay = np.asarray(delays[i])
        this_depth = np.asarray(depths[i])

        this_mean = np.mean(this_delay)
        my_vals = np.where(np.abs(this_delay-this_mean)<cut_vals[i])
        this_depth = this_depth[my_vals]
        this_delay = this_delay[my_vals]
        fitvals = np.polyfit(this_delay,this_depth,1)
        test = np.abs(this_depth-(fitvals[1]+fitvals[0]*this_delay))
        #plt.scatter(this_delay,this_depth,color='darkblue')
        if(linear_fit[i]==0):
            this_delay = this_delay[test<15.0]
            this_depth = this_depth[test<15.0]

        new_delay = []
        new_depth = []
        for j in range(0,len(this_delay)/10):
            new_delay.append(np.mean(this_delay[j*10:j*10+10]))
            new_depth.append(np.mean(this_depth[j*10:j*10+10]))
        new_delay = np.asarray(new_delay)
        new_depth = np.asarray(new_depth)

        #plt.scatter(this_delay,this_depth,color='dodgerblue')
        #plt.plot(new_delay,new_depth,color='darkblue',lw=2.0)
        #plt.show()
        f_int = interpolate.interp1d(new_depth,new_delay)
        temp_delay = np.zeros(len(SPICE_depths))
        for d in range(0,len(SPICE_depths)):
            temp_delay[d]=f_int(SPICE_depths[d])
        interp_delays.update({i:np.asarray(temp_delay)})
    #np.save('all_SPICE_delays_int.npy',interp_delays)
    #np.save('all_SPICE_depths_int.npy',SPICE_depths)


    P = -114.6#-114.78
    Q = -116.6#-117.08
    general_r = 4160.0
    general_z = -194.52
    best_r = 21.11#20.79#22.11
    best_zd = 4.81#4.68#4.74
    best_c = 10.53#12.63#5.0
    all_cables={(8,9):P,(24,25):Q,(8,24):best_c,(9,25):-P+best_c+Q,(8,25):best_c+Q,(9,24):best_c-P}

    for i in pairs:

        good_times = []
        this_r1 = {8:general_r,9:general_r,24:general_r+best_r,25:general_r+best_r}
        this_r2 = {8:general_r,9:general_r,24:general_r+best_r,25:general_r+best_r}
        this_z1 = {8:general_z,9:general_z+29.85,24:general_z+best_zd,25:general_z+best_zd+29.94}
        this_z2 = {8:general_z,9:general_z+29.85,24:general_z+best_zd,25:general_z+best_zd+29.94}

        for zs in SPICE_depths:
            t1 = interp([this_r1[i[0]],zs,this_z1[i[0]]])
            t2 = interp([this_r2[i[1]],zs,this_z2[i[1]]])
            good_times.append(t1-t2-all_cables[i])
        plt.figure(1,facecolor='w')
        this_delay = np.asarray(delays[i])
        this_depth = np.asarray(depths[i])
        plt.scatter(-this_delay,this_depth,color='dodgerblue')
        plt.plot(-1*np.asarray(good_times),SPICE_depths,color='darkblue',lw=3.5)
        #plt.plot(-interp_delays[i],SPICE_depths,lw=2.0)
        plt.xlabel('Time Delay (ns)')
        plt.ylabel('Depth (m)')
        plt.title('Pair:'+str(i))
        plt.show()

    best_chi = 1000000.0


    distances = np.linspace(20.0,22.0,10)
    depth_diffs = np.linspace(3.33,5.0,10)
    cables = np.linspace(5.0,20.0,20)
    cal_distances = np.linspace(25.0,25.0,1)

    chi2 =np.zeros([len(depth_diffs),len(distances),len(cables),len(cal_distances)])
    LoadModel(1.326,1.78,0.0202)
    for rv in range(0,len(distances)):
        print(distances[rv])
        for rc in range(0,len(cal_distances)):

            for z in range(0,len(depth_diffs)):
                for R in range(0,len(cables)):
                    #print(R)
                    all_cables={(8,9):P,(24,25):Q,(8,24):cables[R],(9,25):-P+cables[R]+Q,(8,25):cables[R]+Q,(9,24):cables[R]-P}
                    for i in pairs:
                        this_r1 = {8:general_r,9:general_r,24:general_r+distances[rv],25:general_r+distances[rv]}
                        this_r2 = {8:general_r,9:general_r,24:general_r+distances[rv],25:general_r+distances[rv]}
                        this_z1 = {8:general_z,9:general_z+29.85,24:general_z+depth_diffs[z],25:general_z+depth_diffs[z]+29.94}
                        this_z2 = {8:general_z,9:general_z+29.85,24:general_z+depth_diffs[z],25:general_z+depth_diffs[z]+29.94}

                        counter = 0
                        #print(i,all_cables[i])

                        for zs in SPICE_depths:
                            t1 = interp([this_r1[i[0]],zs,this_z1[i[0]]])
                            t2 = interp([this_r2[i[1]],zs,this_z2[i[1]]])

                            #t1 = trace(this_r1,zs,this_z1)#ch 8 is the closer one
                            #t2 = trace(this_r2,zs,this_z2)#ch 24 is farther
                            #print(rv,z,R,rc)
                            chi2[z,rv,R,rc]=chi2[z,rv,R,rc]+((t1-t2-all_cables[i])[0]-interp_delays[i][counter])**2

                            #print(i,(t1-t2)[0],all_cables[i],interp_delays[i][counter])
                            counter=counter +1
                        """
                        this_r1 = {24:37.88,25:38.5,8:37.88+cal_distances[rc],9:37.88+cal_distances[rc]}
                        this_r2 = {24:37.88,25:38.5,8:37.88+cal_distances[rc],9:37.88+cal_distances[rc]}
                        #print(this_r1[i[0]],this_z1[i[0]],this_r2[i[1]],this_z2[i[1]])
                        t1 = trace(this_r1[i[0]],-173.28,this_z1[i[0]])
                        t2 = trace(this_r2[i[1]],-173.28,this_z2[i[1]])
                        chi2[z,rv,R,rc]=chi2[z,rv,R,rc]+(t1-t2-all_cables[i]-cal_delays[i])**2
                        """
                        #print(i,t1,t2,all_cables[i],cal_delays[i])
                    #print(distances[rv],depth_diffs[z],cables[R],chi2[rv,z,R])
                    if(chi2[z,rv,R,rc]<best_chi):
                        print('')
                        print(distances[rv],depth_diffs[z],cables[R],cal_distances[rc],chi2[z,rv,R,rc])
                        best_chi = chi2[z,rv,R,rc]



def CalculatePropTimes():
    rvals = np.linspace(4150.0,4190.0,80)
    zs_vals = np.linspace(-1450,-1050,20)
    zd_vals = np.linspace(-200.0,-145.0,110)
    prop_times = np.zeros([len(rvals),len(zs_vals),len(zd_vals)])
    LoadModel(1.326,1.78,0.0202)
    for r in range(0,len(rvals)):
        print('')
        print(r)
        for zs in range(0,len(zs_vals)):
            #print(zs)
            for zd in range(0,len(zd_vals)):
                prop_times[r,zs,zd] = trace(rvals[r],zs_vals[zs],zd_vals[zd])#ch 8 is the closer one
    np.save('SPICE_prop_times.npy',prop_times)
    return(prop_times)

def LoadPropTimes():
    prop_times = np.load('SPICE_prop_times.npy')
    rvals = np.linspace(4150.0,4190.0,80)
    zs_vals = np.linspace(-1450,-1050,20)
    zd_vals = np.linspace(-200.0,-145.0,110)
    interpolator = RegularGridInterpolator((rvals,zs_vals,zd_vals),prop_times)
    print(interpolator([4160,-1090.,-180]))
    return(interpolator)
def main():
    #CalculatePropTimes()
    #LoadPropTimes()

    #t0,v0 = LoadSpiceData('5',"3989","0208","2018",int(0),1792,0,1,0,0,0)
    print('')
    """
    t0,v0 = LoadSpiceData('5',"2052","0110","2018",int(0),2048,0,1,1,1,0)
    t1,v1 = LoadSpiceData('5',"2052","0110","2018",int(1),2048,0,1,1,1,0)
    t8,v8 = LoadSpiceData('5',"2052","0110","2018",int(8),2048,0,1,1,1,0)
    t9,v9 = LoadSpiceData('5',"2052","0110","2018",int(9),2048,0,1,1,1,0)
    t16,v16 = LoadSpiceData('5',"2052","0110","2018",int(16),2048,0,1,1,1,0)
    t17,v17 = LoadSpiceData('5',"2052","0110","2018",int(17),2048,0,1,1,1,0)
    t24,v24 = LoadSpiceData('5',"2052","0110","2018",int(24),2048,0,1,0,0,0)
    t24 = t24-t24[0,0]
    t25,v25 = LoadSpiceData('5',"2052","0110","2018",int(25),2048,0,1,1,1,0)

    """
    t0,v0 = LoadSpiceData(int(0),1,1,1)

    t1,v1 = LoadSpiceData(int(1),1,1,1)

    t8,v8 = LoadSpiceData(int(8),1,1,1)
    t9,v9 = LoadSpiceData(int(9),1,1,1)
    t16,v16 = LoadSpiceData(int(16),1,1,1)
    t17,v17 = LoadSpiceData(int(17),1,1,1)
    t24,v24 = LoadSpiceData(int(24),1,0,0)
    #t24 = t24-t24[0,0]
    t24=t24-0.3125 #subtract so that it starts at 0 like the other channels
    t25,v25 = LoadSpiceData(int(25),1,1,1)

    #print(t0,t24)


    delay_dict = {}
    depth_dict = {}
    ts = {0:t0,1:t1,8:t8,9:t9,16:t16,17:t17,24:t24,25:t25}
    vs = {0:v0,1:v1,8:v8,9:v9,16:v16,17:v17,24:v24,25:v25}
    #chs = {0:0,1:1,2:8,3:9,4:16,5:17,6:24,7:25}

    pairs=list(itertools.combinations((0,1,8,9,24,25,16,17),2))

    for i in pairs:
        print(ts[i[0]])
        delays24,depth24 = FindSPICEDelays(ts[i[0]],ts[i[1]],vs[i[0]],vs[i[1]],i[0],i[1])
        delay_dict.update({i:delays24})
        depth_dict.update({i:depth24})
    np.save('delay_dict.npy',delay_dict)
    np.save('depth_dict.npy',depth_dict)
    print('delay dict updated')
    """
    delay_dict_np = np.load('delay_dict.npy',allow_pickle=True)
    depth_dict_np = np.load('depth_dict.npy',allow_pickle=True)
    for i in pairs:
        delay_dict.update({i:delay_dict_np.item().get(i)})
        depth_dict.update({i:depth_dict_np.item().get(i)})
    #print(delay_dict)
    #print(delay_dict.item().get((24,25)))
    #SpiceChi(delay_dict.item().get((24,25)),depth_dict.item().get((24,25)))
    #SpiceChi(delay_dict.item().get((8,9)),depth_dict.item().get((8,9)))
    #SpiceChi(delay_dict[(0,1)],depth_dict[(0,1)])
    #SpiceChi(delay_dict[(16,17)],depth_dict[(16,17)])

    #plt.legend()
    #plt.show()
    TwoStringFitter(delay_dict,depth_dict,pairs,8,24)
    #TwoStringSpiceChi(delay_dict,depth_dict,pairs,8,24)
    print(np.shape(v1))
    ts = {0:t0,1:t1,2:t8,3:t9,4:t16,5:t17,6:t24,7:t25}
    vs = {0:v0,1:v1,2:v8,3:v9,4:v16,5:v17,6:v24,7:v25}
    for i in range(40,1500):
        if(np.max(v0[i])<900):
            plt.figure(1,facecolor='w')
            for j in range(0,8):
                #plt.plot(t0[i]-t0[0,0],v0[i])
                plt.subplot(4,2,j+1)
                plt.plot(ts[j][i],vs[j][i])
            plt.show()
    """

if __name__=="__main__":
   main()
