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
#sys.path.append('../')
from AutomaticLoadData import LoadARACalPulsers
#from TimingCalibration import partial_derivative
#from VoltageCalibration import BlockCorrector, MeanEachBlock
#from AutomaticLoadData import LoadDataFromWeb, DownloadPedFile,LoadFilteredData
import matplotlib
import scipy
import itertools
#from CalPulser import loadvalues
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def AvWaveform(channel):
    print('Averaging waveform for channel', channel)

    num_blocks = 188 #192
    total_samples = 1792#896

    colors =["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

    #time,t_cal,ADC,volt,block_nums =loadvalues(channel,total_samples,num_blocks,rfile,ped_values_ARA)
    #t_cal, volt = LoadDataFromWeb("5337","0529","2019",int(channel),1,1,1,1)
    if(int(channel)==24):
        t_cal,volt=LoadARACalPulsers('3416',int(channel),1,1,0,0)

    elif(int(channel) in [0,1,8,9,16,17,25]):

        t_cal,volt=LoadARACalPulsers('3416',int(channel),1,1,1,1)

    #t_cal = t_cal*1.04
    total_samples = len(t_cal)
    #print(total_samples)

    comp_event = 5
    sig_orig = volt[comp_event]/np.std(volt[comp_event])
    #print(t_cal,sig_orig)
    #plt.plot(t_cal,sig_orig)
    #plt.show()
    #sig_orig = volt_orig_up/(np.std(volt_orig_up))
    averageV = volt[comp_event]
    #print('original length',len(t_cal_up))
    #averageV = RemoveAngularResponse(t_cal_up,averageV,channel,2)
    dt =t_cal[1]-t_cal[0]
    for wf in range(1,len(volt[:,0])):

        sig_comp = volt[wf]/(np.std(volt[wf]))
        #plt.plot(t_cal_up,volt_comp)
        #plt.show()
        #print(len(sig_comp),len(sig_orig))


        cor = 1.0/(float(len(sig_orig)-1))*signal.correlate(sig_orig,sig_comp)
        #print(np.max(cor))


        delay=(np.argmax((cor))-np.size(cor)/2.)*dt
        #print(delay)
        volt_new = np.roll(volt[wf],int(delay/dt))

        averageV = averageV+volt_new

    averageV=averageV/len(volt[:,0])

    return(averageV,t_cal)

def FindDelays(channel,rfile,Av,t_av):
    #ped_values_ARA = np.loadtxt('../run_004804/pedestalValues.run004804.dat')
    #ped_values_ARA = np.loadtxt('../run_005372/pedestalValues.run005372.dat')
    #num_blocks = 188
    total_samples = 2048

    colors =["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
    if(int(channel)==24):
        #t_cal, volt,blocks = LoadDataFromWeb('5',"3416","0602","2018",int(channel),total_samples,1,1,0,0,0)
        #t_cal = t_cal-t_cal[0,0]
        t_cal,volt=LoadARACalPulsers('3416',int(channel),1,1,0,0)
    elif(int(channel) in [0,1,8,9,16,17,25]):
        #t_cal, volt,blocks = LoadDataFromWeb('5',"3416","0602","2018",int(channel),total_samples,1,1,1,1,0)
        t_cal,volt=LoadARACalPulsers('3416',int(channel),1,1,1,1)
    total_samples = len(volt[:,0])

    SigAv = Av/np.std(Av)
    delays = []
    #print(len(volt[:,0]))
    #print(len(t_cal),len(volt[0,:]))
    dt =t_cal[1]-t_cal[0]
    for wf in range(1,len(volt[:,0])):
        #if(block_nums[wf]%2==0):
        #    t_cal_comp=t_cal[:len(volt[wf,:])]
        #else:
        #    t_cal_comp=t_cal[64:64+len(volt[wf,:])]-t_cal[64]

        #f = interp1d(t_cal_comp,volt[wf,:])
        #volt_comp = f(t_cal_even)
        #print(len(t_cal_even),len(volt[wf,:]))
        #t_cal_even = np.arange(t_cal[wf,0],t_cal[wf,-1]-1.0,0.3125)
        #f = Akima1DInterpolator(t_cal[wf,:],volt[wf,:])
        #volt_comp = f(t_cal_even)

        #volt_comp,t_cal_up = signal.resample(volt_comp,len(volt_comp)*10,t=t_cal_even)
        #volt_comp=RemoveAngularResponse(t_cal_up,volt_comp,channel,2)
        sig_comp = volt[wf]/(np.std(volt[wf]))
        #sig_comp = sig_comp/np.max(sig_comp)

        cor = 1.0/(float(len(SigAv)))*signal.correlate(SigAv,sig_comp)
        delay=(np.argmax((cor))-np.size(cor)/2.)*dt+t_cal[0]-t_av[0]
        #print('delay is,',delay)
        delays.append(delay)

        #plt.plot(t_cal,sig_comp)
        #plt.show()


        volt_new = np.roll(volt[wf],int(delay/dt))

    #print(delays)
    #plt.figure(0)
    #plt.hist(delays,bins=50)
    #plt.show()
    #print(len(delays))
    return(delays)


def main():



    Av16,t16 = AvWaveform('16')
    Av24,t24 = AvWaveform('24')
    Av0,t0 = AvWaveform('0')
    Av1,t1 = AvWaveform('1')
    Av8,t8 = AvWaveform('8')
    Av9,t9 = AvWaveform('9')

    Av17,t17 = AvWaveform('17')
    Av25,t25 = AvWaveform('25')




    topv = ['1','9','17','25']
    botv = ['0','8','16','24']
    long_delay = 0.0#132.3
    short_delay = 0.0#18.62

    Avs = {0:Av0,1:Av1,2:Av8,3:Av9,4:Av24,5:Av25,6:Av16,7:Av17}
    ts = {0:t0-long_delay,1:t1-short_delay,2:t8-long_delay,3:t9-short_delay,4:t24-long_delay,5:t25-short_delay,6:t16-long_delay,7:t17-short_delay}

    #Avs = {0:Av0,1:Av1,2:Av8,4:Av24,5:Av25,6:Av16,7:Av17}
    #ts = {0:t0-132.3,1:t1-18.62,2:t8-132.3,4:t24-132.3,5:t25-18.62,6:t16-132.3,7:t17-18.62}
    plt.figure(0,facecolor='w')
    for i in range(0,8):
        print(i)
        #if(i!=4):
        plt.subplot(4,2,i+1)
        plt.plot(ts[i],Avs[i])
        #plt.xlim([-100,200])
        plt.grid()
        #plt.xlim([50,350])
        plt.ylim([-450,450])
        #plt.ylim([-500,500])
        #this_t = ts[i]
        #plt.axvline(this_t[np.argmax(Avs[i])],color='red',lw=3.0)
        #hilbert_env=hilbert(Avs[i])
        #plt.plot(ts[i],np.abs(hilbert_env))
    plt.savefig('plots/AverageCalPulserWFs.pdf')
    plt.show()
    Avs = {0:Av0,1:Av1,8:Av8,9:Av9,24:Av24,25:Av25,16:Av16,17:Av17}

    ts = {0:t0-long_delay,1:t1-short_delay,8:t8-long_delay,9:t9-short_delay,24:t24-long_delay,25:t25-short_delay,16:t16-long_delay,17:t17-short_delay}
    #ts[0]=ts[0]+6.8
    #ts[8]=ts[8]+0.9
    #ts[16]=ts[16]+4.0
    #ts[24]-ts[24]+3.2

    #NOTE: try using zero crossing after biggest peak!
    #Avs = {0:Av0,1:Av1,8:Av8,9:Av9,25:Av25,16:Av16,17:Av17}
    delay0 = np.asarray(FindDelays('0','1402',Av0,t0))
    #plt.figure(0,facecolor='w')
    #plt.hist(delay0,bins=100)
    #plt.show()
    delay1 = np.asarray(FindDelays('1','1411',Av1,t1))
    delay8 = np.asarray(FindDelays('8','1402',Av8,t8))
    delay9 = np.asarray(FindDelays('9','1411',Av9,t9))
    delay16 = np.asarray(FindDelays('16','1402',Av16,t16))
    delay17 = np.asarray(FindDelays('17','1411',Av17,t17))
    delay24 = np.asarray(FindDelays('24','1402',Av24,t24))
    delay25 = np.asarray(FindDelays('25','1411',Av25,t25))


    delays = {0:delay0,1:delay1,8:delay8,9:delay9,24:delay24,25:delay25,16:delay16,17:delay17}

    pairs = list(itertools.combinations(('0','1','8','9','24','25','16','17',), 2))

    #for test in ('0','1','8','9','25','16','17'):
    #    print((1.0/(t0[1]-t0[0]))*1e9)
    #    Avs[int(test)]= butter_lowpass_filter(Avs[int(test)], 400e6, (1.0/(t0[1]-t0[0]))*1e9, order=5)
    print(t0)
    print(t1)
    print(t8)

    #fig,axs = plt.subplots(3,7,sharex=True,facecolor='w')
    plt.figure(1,facecolor='w')
    plt.figure(2,facecolor='w')
    plt.plot(t0,Avs[0]/np.max(Avs[0]))

    cor0 = []
    cor1 = []
    cor2 = []
    counter = 1
    error_dict = {}
    ARA_dict = {}

    for antennas in pairs:
        print('')
        print(ts[int(antennas[0])][0])
        print(ts[int(antennas[1])][0])
        sig0 = Avs[int(antennas[0])]/np.max(Avs[int(antennas[0])])
        sig1 = Avs[int(antennas[1])]/np.max(Avs[int(antennas[1])])
        sig0 = sig0/np.std(sig0)#Avs[int(antennas[0])]/np.std(Avs[int(antennas[0])])
        sig1 = sig1/np.std(sig1)#Avs[int(antennas[1])]/np.std(Avs[int(antennas[1])])
        cor = 1.0/(float(len(Av0)))*signal.correlate(sig0,sig1)
        #print(np.max(cor))
        delay=(np.argmax((cor))-np.size(cor)/2.)*(t0[1]-t0[0])+ts[int(antennas[0])][0]-ts[int(antennas[1])][0]

        diff = delays[int(antennas[0])]-delays[int(antennas[1])]#delay1-delay0
        #diff = diff[diff>(np.mean(diff)-1.0)]
        #diff = diff[diff<(np.mean(diff)+1.0)]
        #print(diff)
        mean,std = norm.fit(diff)
        print('Pairs:',antennas)
        print(delay)

        val = (delay-mean)
        print(val)
        print(np.mean(diff))
        #val = 0
        error_dict.update({(int(antennas[0]),int(antennas[1])):std})
        ARA_dict.update({(int(antennas[0]),int(antennas[1])):val})
        matplotlib.rcParams.update({'font.size': 10})
        binwidth=0.3125/4.0
        plt.figure(1)
        plt.subplot(7,4,counter)
        plt.hist(diff+val,bins=np.arange(min(diff+val), max(diff+val) + binwidth, binwidth),edgecolor='None')
        plt.plot([], [], ' ', label='mean: '+str(round(val,2)))
        plt.plot([], [], ' ', label='std: '+str(round(std,3)))
        plt.legend()

        plt.tick_params(labelsize=10)
        #plt.xlabel(None)
        counter = counter + 1

        if(antennas[0]=='0'):
            plt.figure(2)
            plt.plot(ts[int(antennas[1])]+delay,Avs[int(antennas[1])]/np.max(Avs[int(antennas[1])]),lw=2.0)

    print(ARA_dict)
    print(error_dict)

    np.save('data/ARA_CalPulser_delaydict.npy',ARA_dict)
    plt.figure(1)
    plt.tight_layout()
    plt.savefig('plots/CalPulser_Delays.pdf')
    plt.show()









if __name__=="__main__":
   main()
