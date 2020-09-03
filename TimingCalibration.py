from __future__ import print_function
from ROOT import TCanvas, TGraph
from ROOT import gROOT
import ROOT
import matplotlib.pyplot as plt
from scipy.fftpack import rfft, irfft, fftfreq,fft
from scipy import optimize
from scipy.misc import derivative
import numpy as np
import sys
import math
#sys.path.append('../')
#import global_constants
from math import sin
from pynverse import inversefunc
from AutomaticLoadData import LoadSineWaveData
import matplotlib

matplotlib.rcParams.update({'font.size': 26})
font = {'weight' : 'bold',
        'size'   : 18}
matplotlib.rc('font', **font)

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


def invertedFit(params,tval,vval):
    k = params[0]
    phi = params[1]
    A = params[2]
    T = 1/(k)
    #print(T,v_off,phi,k,A)


    sine = (lambda t: A*np.sin(2*np.pi*k*t-phi))

    t_list = np.linspace(tval-T/2,tval+T/2,100)
    sine_vals = sine(t_list)
    a_max = np.argmax(sine_vals)
    a_min = np.argmin(sine_vals)


    t_max=t_list[a_max]
    t_min=t_list[a_min]


    if(t_min>t_max):
        minval=t_max
        t_max=t_min
        t_min=minval

    #while(t_min>tval):
    #    t_min=tval-0.05
    #while(t_max<tval):
    #    #   print('hello!')
    #    t_max=tval+0.025
    #    #   print(t_max)

    try:
        t_close = inversefunc(sine,y_values=vval,domain=[t_min,t_max])
    except ValueError:
        print(t_max,t_min,vval,tval)
        print(params)

        #plt.figure(1)
        #plt.plot(np.linspace(0,40,150),sine(np.linspace(0,40,150)))
        #plt.show()

    #jitter= t_close-tval
    #if(jitter>

    return(t_close-tval)#used to be t_close-tval


def SineFunc(t,k,phi,A): #time, freq, offset, amplitude
    return A*np.sin(2*np.pi*k*t-phi)

def SineFit(t,v,freq):
    params, params_covariance = optimize.curve_fit(SineFunc,t,v,p0=[freq,np.pi/2.0,350])#,bounds=([-np.inf,-np.inf,200],[np.inf,np.inf,np.inf]))#freq,offset,amplitude,voff
    if(params[2]<0):
        params[2]=np.abs(params[2])
        params[1]=params[1]+np.pi
    params[1]=params[1]%(np.pi*2)
    while(params[1]<0):
        params[1]=params[1]+np.pi*2.0
    return(params)


def HistPlotter2D(sample,jitter):
    sample_even=[]
    jitter_even = []
    sample_odd=[]
    jitter_odd = []

    sample_1d = 4
    jitter_1d = []

    for j in range(0,len(sample)):
        if(sample[j]==sample_1d):
            jitter_1d.append(jitter[j])

        if(sample[j]%2==0):#if even
            sample_even.append(sample[j])
            jitter_even.append(jitter[j])
        else:
            sample_odd.append(sample[j])
            jitter_odd.append(jitter[j])
    #print(sample_even)
    #print(jitter_even)
    plt.figure(4,facecolor='w')
    plt.hist2d(sample_even,jitter_even,bins=(128,128),cmap=plt.cm.jet,range=np.array([(0.0,128.0),(-1.0,1.0)]))
    plt.title('Even Samples')

    plt.figure(5,facecolor='w')
    plt.hist2d(sample_odd,jitter_odd,bins=(128,128),cmap=plt.cm.jet,range=np.array([(0.0,128.0),(-1.0,1.0)]))
    plt.title('Odd Samples')

    plt.figure(3)
    plt.hist(jitter_1d,bins=25)
    plt.xlabel('Jitter (ns)')
    plt.xlim([-1.0,1.0])


    plt.show()
    return()

def AddOffsets2(t,v,freq,odds):
    #Goals:
    #Fit even blocks to sine waves
    #Fit odd blocks to sine waves
    #compare phase differences between two transitions
    print(odds)
    e_freq = []
    o_freq = []
    #loop over all events:
    for j in range(0,len(v[:,0])):
        blocks_per_event = int(len(v[0,:])/64)
        #print('total blocks:', blocks_per_event)
        for i in range(0,blocks_per_event):
            #print(i*64,i*64+64,15*64)
            params = SineFit(t[1:64:2],v[j,i*64:i*64+64:2],freq)

            if i%2==0:#even block
                e_freq.append(params[0])
            else:
                o_freq.append(params[0])
            #print(params)
    histogram([e_freq,o_freq],'Frequency')


    t_scaled_e1 = t[:64]*np.mean(e_freq)/freq
    t_scaled_o = t[64:128]*np.mean(o_freq)/freq
    t_scaled_e2 = t[128:192]*np.mean(e_freq)/freq

def AddOffsets(t,v,freq,odds,old_mean_o2e,old_mean_e2o):
    print(len(v[0,:]))

    e_freq = []
    o_freq = []
    for j in range(0,len(v[:,0])):
        val=0
        for i in range(0,14):
            params = SineFit(t[odds[:32]],v[j,odds[val:val+32]],freq)
            if i%2==0:#even block
                e_freq.append(params[0])
            else:
                o_freq.append(params[0])
            val = val+32

    t_scaled_e1 = t[:64]*np.mean(e_freq)/freq
    t_scaled_o = t[64:128]*np.mean(o_freq)/freq
    t_scaled_e2 = t[128:192]*np.mean(e_freq)/freq

    e2o_diff = []
    o2e_diff = []
    #plt.figure(0)
    #plt.plot(even_time[odds[:32]],v[0,odds[0:32]])
    #plt.show()

    for j in range(0,len(v[:,0])):
        val=0
        for i in range(0,13):
            #print(i)
            if(i%2==0):
                params = SineFit(t_scaled_e1[odds[:32]],v[j,odds[val:val+32]],freq)
                #print(params[0])
                e_offset=(params[1])
                params = SineFit(t_scaled_o[odds[:32]],v[j,odds[val+32:val+64]],freq)
                o_offset=(params[1])
                diff = e_offset-o_offset
                while(diff>np.pi):
                    diff=diff-2*np.pi
                while(diff<-1*np.pi):
                    diff=diff+2*np.pi
                #print(e_offset,o_offset)
                e2o_diff.append(((diff))/(2*np.pi*freq))
            else:
                #print(oe_time[:32])
                #print(len(oe_time[:32]))
                params = SineFit(t_scaled_o[odds[:32]],v[j,odds[val:val+32]],freq)
                o_offset=(params[1])
                params = SineFit(t_scaled_e2[odds[:32]],v[j,odds[val+32:val+64]],freq)
                e_offset=(params[1])
                diff = e_offset-o_offset
                while(diff>np.pi):
                    diff=diff-2*np.pi
                while(diff<-1*np.pi):
                    diff=diff+2*np.pi
                o2e_diff.append(((diff))/(2*np.pi*freq))
            val = val+32

    #o2e_diff = reject_outliers(np.asarray(o2e_diff))
    #histogram([e2o_diff,o2e_diff],'Wrap Around Time (ns)')
    #Scale time so that offset is taken into account
    e2o_mean = np.abs(np.mean(e2o_diff)+old_mean_e2o)
    o2e_mean = np.abs(np.mean(o2e_diff)+old_mean_o2e)
    print(e2o_mean,o2e_mean,old_mean_e2o,old_mean_o2e)
    if(e2o_mean>1.0):
        e2o_mean=0.31
    if(o2e_mean>1.0):
        o2e_mean=0.31
    print('the offsets are:', e2o_mean, o2e_mean)

    #histogram([e2o_diff,o2e_diff],'time offset')
    t_updated = np.zeros(896)
    val = 0
    for i in range(0,7):
        if i==0:
            t_updated[val:val+64]=t_scaled_e1
            t_updated[val+64:val+128]=t_scaled_o+e2o_mean
        else:
            t_updated[val:val+64]=t_scaled_e1+o2e_mean+t_updated[val-1]
            t_updated[val+64:val+128]=t_scaled_o+e2o_mean+t_updated[val-1]
        #else:
        #    t_updated[val:val+64]=even_time+(o2e_mean+t_updated[val-1])
        #t_updated[val+64:val+128]=odd_time+(e2o_mean+t_updated[val+64-1])
        val = val+128
    print(t_updated)
    return(t_updated,o2e_mean,e2o_mean)

    #print(t)
def histogram(vals,string):
    fig, axs = plt.subplots(2,1,facecolor='w')
    counter=0
    for ax in axs.reshape(-1):
        print('success')
        ax.hist(vals[counter],color='navy',edgecolor='none',bins=20)
        ax.axvline(x=np.mean(vals[counter]),color='red',ls='-',linewidth=2.0)
        #ax.text(230,250,"mean (MHz):  "+str(round(np.mean(freq_array[counter]*1000),2)))
        #ax.set_xlim(200,250)
        #ax.set_ylim(0,300)
        ax.set_xlabel(string)
        ax.set_ylabel('Counts')
        counter = counter +1
    plt.show()

def SinePlotter(t,v,params,sample):
    plt.figure(10)
    plt.scatter(t,v[sample,:],label='Data')
    #print(t)
    #print(v[sample,:])
    t_up = np.linspace(t[0],t[-1],len(t)*50)
    plt.plot(t_up,SineFunc(t_up,params[sample,0],params[sample,1],params[sample,2]),label='Best Fit: '+str(np.round(params[sample,0]*1000,3))+' MHz')
    plt.grid()
    #plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('ADC')
    plt.show()

def SlopeFinder(t,v,sample):
    if sample==0:
        slope = (v[sample+2]-v[sample])/(t[sample+2]-t[sample])
    else:
        slope=(v[sample]-v[sample-2])/(t[sample]-t[sample-2])

    return slope

def CorrectTimingSample(rootfile,channel,freq,t_cal,station):
    wf_len = 896
    
    pedFile = '/home/kahughes/PA_Analysis/data/pedFiles/pedFile_'+str(rootfile)+'.dat'

    print(pedFile)

    all_times,volt,blocks = LoadSineWaveData(station,rootfile,pedFile,channel,kPed=1,kTime=0,kVolt=0)

    time = all_times[0]-all_times[0][0] #uncalibrated time

    print('number of events is', np.shape(volt))
    num_blocks=len(volt[:,0])

    best_params = np.zeros([num_blocks,4])

    odds = np.linspace(1,wf_len-1,int(wf_len/2),dtype=int)
    evens = np.linspace(0,wf_len-2,int(wf_len/2),dtype=int)

    odd_params=np.zeros([num_blocks,3])

    t_cal = np.zeros(128)


    if(t_cal[5]==0.0):
        print('clearing out old t_cal')
        t_cal=time[:128]
        print(t_cal)

    #set calibrated time array to be equal to uncalibrated time array, so length is the same

    t_cal_full = time 
    #print('t_cal before is', t_cal_full)
    #print(np.shape(t_cal_full),np.shape(volt))
    odd_mean= 0.0
    even_mean= 0.0
    for l in range(0,4):
        print('loop number ', l)
        #print(np.shape(t_cal_full),np.shape(volt))
        #First fix offsets between blocks, as that can be larger than the channel to channel fixes.
        #if l==0:
        #    t_cal_full,odd_mean,even_mean=AddOffsets(t_cal_full,volt,freq,odds,odd_mean,even_mean)
            #t_cal_full,odd_mean,even_mean=AddOffsets2(t_cal_full,volt,freq,odds)
        #Fit each waveform to a sine wave
        volt=volt[:,:len(t_cal_full)]
        #print(np.shape(t_cal_full),np.shape(volt))


        #STEP ONE: For each event, calculate the best sine wave fit. Then remove outlier frequencies and calculate the mean frequency.
        #Since we know the expected frequency from the lab data, we can then correct the overall time step.

        for i in range(0,num_blocks):
            odd_params[i,:]=SineFit(t_cal_full[odds],volt[i,odds],freq)
            #SinePlotter(t_cal_full,volt,odd_params,i)

        #plt.scatter(t_cal_full[odds],volt[i,odds])
        #t_up = np.linspace(t_cal_full[0],t_cal_full[-1],500)
        #plt.plot(t_up,SineFunc(t_up,odd_params[i,0],odd_params[i,1],odd_params[i,2]))
        #plt.show()

        """
        if(l>-1):
            plt.hist(odd_params[:,0]*1000,bins=25)
            #plt.grid()
            plt.axvline(np.mean(odd_params[:,0])*1000,color='red',label='Mean')
            plt.axvline(218.0, color='black',label='Input Frequency')
            plt.xlabel('Best Fit Frequency (MHz)')
            plt.legend()
            plt.show()
        
            SinePlotter(t_cal_full,volt,odd_params,4)
        """
        freq_no_outliers = reject_outliers(np.asarray(odd_params[:,0]))
        mean_freq = np.mean(freq_no_outliers)
        print('mean frequency is', mean_freq)
        #histogram([odd_params[:,1],even_params[:,1]],'')

        #Scale timing to reflect true frequency
        t_cal_full=t_cal_full*mean_freq/freq

        print('Fitting to sine:')
        #Re-fit using new time
        for i in range(0,num_blocks):

            odd_params[i,:]=SineFit(t_cal_full[odds],volt[i,odds],freq)
            #SinePlotter(t_cal_full,volt,odd_params,i)
            #print('here')
            #print(l)
            #if(l>0):
                #print("here!")
                #SinePlotter(t_cal_full,volt,odd_params,i)
                #plt.show()

        t_cal = t_cal_full[:128]

        jitter_array = []
        sample_array = []
        slope_array = []
        jitter_slope = []
        new_spacing = np.zeros(128) #spacing between 0 and 1, 1 and 2, etc.
        #Go through each sample and correct based on fit

        if(int(channel)==24):
            cutval = 5.0
        else:
            cutval = 30.0

        print('here is the slow part')
        for k in range(0,896):
            counter = 0
            for i in range(0,num_blocks):
                if(np.abs(volt[i,k])<cutval and (freq-odd_params[i,0])<0.002):# and np.abs(odd_params[i,2]>200)):
                    try:
                        invert_fit = invertedFit(odd_params[i,:],t_cal_full[k],volt[i,k])
                        jitter_array.append(invert_fit)
                        sample_array.append(k%128)
                        counter = counter+1
                    except:
                        print('error in finding inverse!')

            t_cal_full[k]=t_cal_full[k]+np.mean(jitter_array[-counter:])

            if(k>0):
                new_spacing[k%128]=new_spacing[k%128]+t_cal_full[k]-t_cal_full[k-1]


        new_spacing[1:]=new_spacing[1:]/7.0
        new_spacing[0]=new_spacing[0]/6.0
        #print('spacing is', new_spacing)



        for i in range(0,896):
            if(i==0):
                t_cal_full[i]=0.0
            else:
                t_cal_full[i]=t_cal_full[i-1]+new_spacing[(i)%128]


        #print('final t_cal is',t_cal_full)

        t_cal=t_cal_full[:128]

        #if(l>0):
        #    SinePlotter(time[odds],volt[:,odds],odd_params,5)


        """
        plt.figure(6,facecolor='w')
        plt.hist2d(slope_array,jitter_slope,bins=(250,128),cmap=plt.cm.jet,range=np.array([(-1100.0,1100.0),(-1.0,1.0)]))
        plt.title('Even Samples vs Slope')
        plt.show()
        """


        if(l<1):
            np.save('data/ARA'+str(station)+'_cal_files/samples_'+rootfile+'_'+channel+'first.npy',np.asarray(sample_array))
            np.save('data/ARA'+str(station)+'_cal_files/jitter_'+rootfile+'_'+channel+'first.npy',np.asarray(jitter_array))
        #if(l>-1):
    HistPlotter2D(sample_array,jitter_array)
    #print('final t_cal is', t_cal_full)
    np.save('data/ARA'+str(station)+'_cal_files/t_cal_'+channel+'.npy',t_cal_full)
    np.save('data/ARA'+str(station)+'_cal_files/samples_'+channel+'final.npy',np.asarray(sample_array))
    np.save('data/ARA'+str(station)+'_cal_files/jitter_'+channel+'final.npy',np.asarray(jitter_array))
    #HistPlotter2D(sample_array,jitter_array)
    #print('t_cal is', t_cal)
    return(t_cal)

def main():


    #rootfile='1402'
    channel = str(sys.argv[1])#'0'
    station = str(sys.argv[2])


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
    #rootfiles = ['1402']
    #rootfile = '1404'
    #if(int(channel)==24):
    #    rootfiles = ['1411']



    #sample_final=np.load('ARA'+str(station)+'_cal_files/samples_'+rootfile+'_'+channel+'final.npy')
    #jitter_final = np.load('ARA'+str(station)+'_cal_files/jitter_'+rootfile+'_'+channel+'final.npy')

    #sample_first=np.load('ARA'+str(station)+'_cal_files/samples_'+rootfile+'_'+channel+'first.npy')
    #jitter_first = np.load('ARA'+str(station)+'_cal_files/jitter_'+rootfile+'_'+channel+'first.npy')
    #HistPlotter2D(sample_first,jitter_first)
    #HistPlotter2D(sample_final,jitter_final)


    freqs = [0.218,0.353,0.521,0.702]

    cal_t= np.zeros(128)

    #jitter = np.load('jitter_1404_3final.npy')
    #samples = np.load('samples_1404_3final.npy')
    #HistPlotter2D(samples,jitter)
    #CorrectTimingSample(rootfiles[0],channel,freqs[0],cal_t,station)

    for a in range(0,1):
        #exists = os.path.isfile('ARA'+str(station)+'_cal_files/jitter_'+rootfiles[a]+'_'+channel+'final.npy')
        #if(exists):
        #    print('file exists!')
        #else:
        #try:
        CorrectTimingSample(rootfiles[a],channel,freqs[a],cal_t,station)
        #except:
        #    print('Error')

    #average all results together
    #average_tcals('ARA'+str(station)+'_cal_files/',channel,rootfiles)

    #account for wrap around time
    #FindWrapAround(rootfiles[1],channel,freqs[1])

if __name__=="__main__":
   main()
