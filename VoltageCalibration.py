from __future__ import print_function
from ROOT import TCanvas, TGraph
from ROOT import gROOT
import ROOT
import os
import matplotlib.pyplot as plt
from scipy import signal,stats
from scipy.fftpack import rfft, irfft, fftfreq
from scipy.signal import iirfilter,lfilter,butter
from scipy import signal
from scipy.fftpack import fft
from scipy import optimize
from scipy.misc import derivative
import numpy as np
import sys
import math
from math import sin
from array import array
from pynverse import inversefunc
import matplotlib
from AutomaticLoadData import LoadSineWaveData


font = {'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)
matplotlib.rcParams['figure.figsize'] = [12, 8]



def SineFunc(t,k,phi,A): #time, freq, offset, amplitude
    return A*np.sin(2.0*np.pi*k*t-phi)


def SineFit(t,v,freq,A):

    test = lambda t,k,phi: SineFunc(t,k,phi,A)
    #print(test(3,5,6))

    params, params_covariance = optimize.curve_fit(lambda t,k,phi: SineFunc(t,k,phi,A),t,v,p0=[freq,np.pi/2.0])#,bounds=([0.216,-np.inf,200],[0.220,np.inf,np.inf]))#freq,offset,amplitude,voff


    #if(params[2]<0):
        #params[2]=np.abs(params[2])
    #   params[1]=params[1]+np.pi
    params[1]=params[1]%(np.pi*2)
    while(params[1]<0):
        params[1]=params[1]+np.pi*2.0
    params = np.append(params,[A])
    #print(params)
    return(params)

def Cubic(a,p0,p1,p2,p3):
    return(p0*a**3+p1*a**2+p2*a+p3)

def plot_voltage(t,adc,p_pos,p_neg,block_num,freq,A,channel):


    start_block = int(block_num)
    #print('starting block is',start_block)
    length = len(adc)
    #print(len(adc))
    piece = np.linspace(start_block*64,(start_block*64+length)%32768,length,dtype=int)
    v = np.zeros(length)

    for k in range(0,length):

        if(k%64==0 and k>0):
            #print('here!')
            start_block=(start_block+1)%512
        #print(start_block,start_block*64+i%64)
        if(adc[k]>0):
            v[k]=Cubic(adc[k],p_pos[start_block*64+k%64,0],p_pos[start_block*64+k%64,1],p_pos[start_block*64+k%64,2],p_pos[start_block*64+k%64,3])
        if(adc[k]<0.0):
            v[k]=Cubic(adc[k],p_neg[start_block*64+k%64,0],p_neg[start_block*64+k%64,1],p_neg[start_block*64+k%64,2],p_neg[start_block*64+k%64,3])



    #vp =Cubic(adc[adc>0],p_pos[piece,0],p_pos[piece,1],p_pos[piece,2],p_pos[piece,3])
    #vn =Cubic(adc[adc<0],p_neg[piece,0],p_neg[piece,1],p_neg[piece,2],p_neg[piece,3])
    #v = np.extend(vp,vn)
    #print(v)
    params = SineFit(t,v,freq,A)

    plt.figure(0,facecolor='w')
    plt.scatter(t[1::2],adc[1::2],color='dodgerblue',)
    #plt.scatter(t[1::2],v[1::2],color='maroon')
    plt.scatter(t[1::2],v[1::2],color='maroon')
    plt.xlabel('Time (ns)')
    plt.ylabel('Voltage (mV)')
    #plt.plot(t,v,color='maroon')
    t_up = np.linspace(0,300.0,2000)
    plt.plot(t_up,SineFunc(t_up,params[0],params[1],params[2]),color='maroon')
    plt.grid()
    plt.legend()
    #print(params[0],params[1],params[2])
    plt.tight_layout()
    plt.savefig('plots/voltage_plot_Chan_'+str(channel)+'_block'+str(int(block_num))+'.pdf')
    #plt.show()
    plt.clf()

def plot_cubic(this_ADC,this_volt,mean_ADC,mean_volt,p_pos,p_neg,i,meanval,channel):
    colors =["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]

    deg = len(p_pos[i,:])-1

    x = np.sort(this_ADC[this_ADC>=-meanval])
    #x=np.append(x,x_int[0])
    x = np.sort(x)
    xn = np.sort(this_ADC[this_ADC<=-meanval])
    #xn=np.append(xn,x_int[0])

    deg = len(p_pos[i,:])
    print('degree is',deg)
    p_p = np.zeros(len(x))
    p_n = np.zeros(len(xn))
    for p in range(0,deg):
        print(p_pos[i,deg-p-1])
        print(deg-p)
        p_p = p_p+p_pos[i,deg-p-1]*x**p
        p_n = p_n+p_neg[i,deg-p-1]*xn**p


    plt.figure(1,facecolor='w')
    plt.scatter(this_ADC,this_volt,s=2,color='black')
    #plt.scatter(mean_ADC,mean_volt,s=10,color='red')
    plt.plot(x,p_p,color=colors[0],lw=2.0)
    plt.plot(xn,p_n,color=colors[1],lw=2.0)
    plt.xlim([-600,600])
    #plt.plot(x,np.multiply(p_pos[i,0],x**3)+np.multiply(p_pos[i,1],x**2)+np.multiply(p_pos[i,2],x)+p_pos[i,3],color=colors[0],lw=2.0)
    #plt.plot(xn,np.multiply(p_neg[i,0],xn**3)+np.multiply(p_neg[i,1],xn**2)+np.multiply(p_neg[i,2],xn)+p_neg[i,3],color=colors[1],lw=2.0)
    plt.xlabel('ADC Counts')
    plt.ylabel('Voltage (mV)')
    plt.grid()
    plt.tight_layout()
    plt.savefig('plots/sample_ADC_v_'+str(channel)+'_'+str(i)+'.pdf')
    #plt.show()
    plt.clf()


def CorrectVoltage(station,files, channel,freq):

    #pedFile = 'pedFile_'+str(files)+'.dat'
    pedFile = 'pedFiles/pedFile_'+str(files)+'.dat'


    ADC_list = [[] for i in range(32768)]#each entry = block_number*64 +sample%64
    volts_list = [[] for i in range(32768)]
    colors =["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]



    print('Loading RootFiles')
    avg_ADC = []

    #choose amplitude to use:
    A = 445.0



    #load calibrated times
    all_times,ADC,blocks = LoadSineWaveData(station,files,pedFile,int(channel),kPed=1,kTime=1,kVolt=0)

    #can also load uncalibrated times if you want to make plots comparing cal vs. uncal
    times,ADC_raw,block_nums = LoadSineWaveData(station,files,pedFile,int(channel),kPed=1,kTime=0,kVolt=0)



    total_events = len(all_times[:,0])
    print('number of events:',total_events)

    total_samples = len(all_times[0,:])

    odds = np.linspace(1,total_samples-1,int(total_samples/2),dtype=int)
    evens = np.linspace(0,total_samples-2,int(total_samples/2),dtype=int)

    #array to hold fit parameters for each event
    odd_params=np.zeros([total_events,3])


    #Loop through all events and fit data to a sine wave:
    for i in range(0,total_events):
        if(i%100==0):
            print(i)

        odd_params[i,:] = SineFit(all_times[i,odds],ADC[i,odds],freq,A)

    t_cal = all_times[0]


    print(block_nums)
    print(np.shape(times))
    print(np.shape(ADC))


    plot_block = []
    plot_adc = []

    #now loop through all events, determine what volt value is expected based on the fit, 
    #and save both that value and the original ADC value in the array associated with that sample.

    for i in range(0,total_events):
        my_block=int(block_nums[i]) #what is the first block of the event
        #print(my_block)
        for j in range(0,total_samples):#896
            if(j%64==0 and j>0):
                #print(j,my_block)
                my_block=(my_block+1)%512 #if we're on a new block, move to the next block (and restart at 0 after block 511)

            volt_val = SineFunc(t_cal[j],odd_params[i,0],odd_params[i,1],A) #determine was the volt value should be


            if( (ADC[i,j]*volt_val>-5000.0 )): #to catch huge outliers (mainly for HPOL channels):
                ADC_list[(my_block*64+j%64)].append(ADC[i,j])#add original adc value to array
                volts_list[(my_block*64+j%64)].append(volt_val) #add volt value to array
                plot_block.append(my_block*64+j%64)
                plot_adc.append(ADC[i,j])

    #np.save('ADClist.npy',ADC_list)
    #np.save('voltslist.npy',volts_list)


    #ADC_list= np.load('ADClist.npy')
    #volts_list=np.load('voltslist.npy')

    tot = 32768
    counts = 0
    degree = 3
    p_pos=np.zeros([32768,degree+1])
    p_neg=np.zeros([32768,degree+1])
    chi2_p = np.zeros([32768])
    chi2_n = np.zeros([32768])

    zero_vals = np.zeros([32768])

    extra_peds = []

    C1 = [0,1,8,9,16,17,24,25]
    C2 = [3,2,11,10,19,18,27,26]

    if(float(channel) in C1):
        max_range = 32768
    else:
        max_range = 16384


    for i in range(0,max_range):
        if(float(channel) in C1):
            i = i
        else:
            i=i*2+1

        if(i%1000==1):
            print(i)

        this_ADC = np.asarray(ADC_list[i])
        #print('average ADC is', np.mean(this_ADC))
        ADC_mean = np.mean(this_ADC)
        this_volt = np.asarray(volts_list[i])


        #Find x-intercept by finding a linear fit near the origin and solving for the x intercept:
        lin_ADC = []
        lin_volt = []
        for j in range(0,len(this_ADC)):
            if (np.abs(this_ADC[j])<200 and np.abs(this_volt[j])<200):
                lin_ADC.append(this_ADC[j])
                lin_volt.append(this_volt[j])

        myslope,intercept,r,p,stderr = stats.linregress(np.asarray(lin_ADC),np.asarray(lin_volt))

        if(myslope<.1): #If the fit doesn't look linear in a good way, just set origin as interception point.
            zero_vals[i]=0
            x_inter = 0
            intercept = 0
        else: #But, if linear fit looks good, correct the ADC values so that the intercept goes through 0
            zero_vals[i] = intercept/myslope
            x_inter = -intercept/myslope
            this_ADC = this_ADC-x_inter
            intercept=0.0


        if(channel in C2):
            p_pos[i-1,:]=p_pos[i,:]
            p_neg[i-1,:]=p_neg[i,:]



        meas_volt = []
        mean_ADC = []
        var_volt = []
        tmin = np.min(this_ADC)
        tmax = tmin+5.0


        #sort ADC and volts in order:
        sorted = np.argsort(this_ADC)
        sort_ADC = this_ADC[sorted]
        sort_volt = this_volt[sorted]

        #calculate average ADC and volt values for steps of size "my_spacing". This is to manage the loop-like behavior.
        my_spacing = 5
        for k in range(0,int(len(this_ADC)/my_spacing)):
            meas_volt.append(np.mean(sort_volt[k*my_spacing:k*my_spacing+my_spacing]))
            mean_ADC.append(np.mean(sort_ADC[k*my_spacing:k*my_spacing+my_spacing]))
            var_volt.append(np.var(sort_volt[k*my_spacing:k*my_spacing+my_spacing]))


        mean_ADC = np.asarray(mean_ADC)
        meas_volt = np.asarray(meas_volt)
        var_volt = np.asarray(var_volt)


        #add fake datapoints at (0,0) so that both positive and negative solutions go through the same point.
        mean_ADC=np.append(mean_ADC,np.zeros(100))
        meas_volt=np.append(meas_volt,np.zeros(100)+intercept)

        try: #fit cubic for positive and negative half of data
            p_pos[i,:] = np.polyfit(mean_ADC[mean_ADC>=0],meas_volt[mean_ADC>=0],degree)
            p_neg[i,:] = np.polyfit(mean_ADC[mean_ADC<=0],meas_volt[mean_ADC<=0],degree)
        except: #ValueError, RankError: basically if there's a problem, set this fit to be the same as the neighboring odd/even value.
            p_pos[i,:] = p_pos[i-2,:]
            p_neg[i,:] = p_neg[i-2,:]


        mean_ADC = mean_ADC[:-100]
        meas_volt = meas_volt[:-100]

        mean_ADC = mean_ADC[~np.isnan(var_volt)]
        meas_volt = meas_volt[~np.isnan(var_volt)]
        var_volt= var_volt[~np.isnan(var_volt)]

        mean_ADC = mean_ADC[np.nonzero(var_volt)]
        meas_volt= meas_volt[np.nonzero(var_volt)]
        var_volt = var_volt[np.nonzero(var_volt)]

        #calculate chi2 value for both fits:
        pred_v = Cubic(mean_ADC[mean_ADC<=0],p_neg[i,0],p_neg[i,1],p_neg[i,2],p_neg[i,3])
        chi2_n[i]=np.sum(np.divide((meas_volt[mean_ADC<=0]-pred_v)**2,var_volt[mean_ADC<=0]))/(len(var_volt[mean_ADC<=0]))
        pred_v = Cubic(mean_ADC[mean_ADC>=0],p_pos[i,0],p_pos[i,1],p_pos[i,2],p_pos[i,3])
        chi2_p[i]=np.sum(np.divide((meas_volt[mean_ADC>=0]-pred_v)**2,var_volt[mean_ADC>=0]))/(len(var_volt[mean_ADC>=0]))


        #This next part of the code will make plots for selected samples, plot a distribution of chi2 values so far, and plot an example corrected waveform.

        
        if((i%960==0 and i>0)):
            print(i)
            print('chi2 is ', chi2_p[i])
            plot_chip = chi2_p[:i]
            plot_chin = chi2_n[:i]
            plt.figure(0)
            plt.hist(plot_chip[plot_chip<100],bins=100,alpha=0.5,label="Positive Chi2")
            plt.hist(plot_chin[plot_chin<100],bins=100,alpha=0.5,label='Negative Chi2')
            plt.xlabel('Chi2 Value')
            plt.ylabel('Counts')
            plt.legend()
            plt.title('Chi2 Distributions')
            plt.tight_layout()
            plt.savefig('plots/Chi2_distributions_'+str(channel)+'.pdf')
            plt.show()
            plt.clf()

            plot_cubic(this_ADC,this_volt,mean_ADC,meas_volt,p_pos,p_neg,i,0.0,channel)
            #print('tcal is ', t_cal,t_cal[1]-t_cal[0],len(t_cal))
            my_ind = np.where(block_nums==counts)
            #print('this index is', my_ind)
            if(my_ind[0][0]!=None):
                plot_voltage(t_cal,ADC[my_ind[0][0],:],p_pos,p_neg,block_nums[my_ind[0][0]],freq,A,channel)
            counts = counts +14
        



    np.save('data/ARA'+str(station)+'_cal_files/p_pos_'+channel+'.npy',p_pos)
    np.save('data/ARA'+str(station)+'_cal_files/p_neg_'+channel+'.npy',p_neg)
    np.save('data/ARA'+str(station)+'_cal_files/chi2_pos_'+channel+'.npy',chi2_p)
    np.save('data/ARA'+str(station)+'_cal_files/chi2_neg_'+channel+'.npy',chi2_n)
    np.save('data/ARA'+str(station)+'_cal_files/zerovals_'+channel+'.npy',zero_vals)


def main():

    channel = str(sys.argv[1])#'0'
    station = str(sys.argv[2])
    freqs = [0.218,0.353,0.521,0.702]


    N1 = [0,3,8,11,16,19,24,27]
    N2 = [1,2,9,10,17,18,25,26]

    if(int(channel) in N1):
        rootfiles = ['1402', '1403','1404','1405']
    if(int(channel) in N2):
        rootfiles = ['1411','1412','1413','1414']


    N1 = [0,3,8,11,16,19,24,27]
    N2 = [1,2,9,10,17,18,25,26]
    N_special = [9,16,24,25]
    if(station=='5'):
        if(int(channel) in N1):
            rootfile='1402'
            rootfiles = ['1402', '1403','1404','1405']
        if(int(channel) in N2):
            rootfile='1411'
            rootfiles = ['1411','1412','1413','1414']
    if(station=='4'):
        if(int(channel) in N1 and int(channel) not in N_special):
            rootfile='2829'
            rootfiles = ['2829', '2830','2831','2832']
        if(int(channel) in N2 and int(channel) not in N_special):
            rootfile='2840'
            rootfiles = ['2840','2841','2842','2843']
        if(int(channel)in N_special):
            rootfiles = ['2855','2856']


    #Currently, only use Frequency = 218 MHz, but could use other frequencies if desired.
    for a in range(0,1):
        CorrectVoltage(station,rootfiles[a],channel,freqs[a])



if __name__=="__main__":
   main()
