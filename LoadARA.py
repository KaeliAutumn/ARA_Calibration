from __future__ import print_function
import os
from os import path
import ROOT
import ctypes
import matplotlib.pyplot as plt
from scipy import signal,stats
from scipy.signal import hilbert,butter,lfilter
import scipy.optimize as opt
from scipy.interpolate import interp1d, Akima1DInterpolator
import numpy as np
import sys
import math
import matplotlib
import itertools
sys.path.append('../')
import global_constants
f#rom load import nuphase_data_reader
f#rom load.Calibrate import Calibrator,LoadNewPed,RemoveFirstBlock,RemoveBackwardsSamples
#import BeaconTau as bt
#import beacon_data_reader
import requests
from bs4 import BeautifulSoup
from scipy.fftpack import fft
import datetime as dt
from scipy.spatial.distance import pdist

import time
ROOT.gSystem.Load(global_constants.ARASIM_LOC+"libAra.so")



def LoadARACalPulsers(run,channel,kCalPulser,kPed,kTime,kVolt):
	#t_cal,volt=LoadFilteredData('5','3416','0602','2018',int(channel),total_samples,1,1,0,0)
	length = 1792
	station= '5'
	#run = '3416'
	date_list=np.load(global_constants.DIR_LOC+'src/load/date_list.npy',allow_pickle='TRUE').item()
	#print(date_list)
	t_file = global_constants.DIR_LOC+"data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
	v_file = global_constants.DIR_LOC+"data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
	#b_file = "SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy"
	if(os.path.isfile(t_file) and os.path.isfile(v_file)):
		t=np.load(t_file)
		v=np.load(v_file)
		#b = np.load(b_file)
		return(t,v)
	else:
		ped_values,address_ped = LoadNewPed(run)
		#test = ROOT.TWebFile.Open("http://icecube:skua@sundog.uchicago.edu/convey/data/exp/ARA/"+year+"/filtered/L0/ARA05/"+date+"/run"+run+"/event"+run+".root")
		test = ROOT.TWebFile.Open(date_list[int(run)])
		calibrator = ROOT.AraEventCalibrator.Instance()
		eventTree = test.Get("eventTree")

		rawEvent = ROOT.RawAtriStationEvent()
		eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))
		totalEvents = eventTree.GetEntries()
		print('total events:', totalEvents)


		t_univ = np.arange(21.0,length*0.3125-1.0,0.03125)

		all_volts = np.zeros([totalEvents,len(t_univ)])
		all_t=np.zeros([totalEvents,len(t_univ)])
		all_blocks=np.zeros([totalEvents])+701

		print('here is t_univ:', t_univ)
		for i in range(0,totalEvents):#totalEvents):

		    eventTree.GetEntry(i)
		    if(rawEvent.isCalpulserEvent()==0 and kCalPulser==1): #if not a cal pulser and we want cal pulsers, go to next event
		        continue
		    if(rawEvent.isCalpulserEvent()==1 and kCalPulser==0):
		        continue
		    #print(i)
		    usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
		    gr1 = usefulEvent.getGraphFromElecChan(channel)

		    t_buff = gr1.GetX()
		    v_buff = gr1.GetY()
		    n = gr1.GetN()
		    t_buff.SetSize(n)
		    v_buff.SetSize(n)
		    t = np.frombuffer(t_buff,dtype=float,count=-1)
		    v = np.frombuffer(v_buff,dtype=float,count=-1)


		    block_number = rawEvent.blockVec[0].getBlock()

		    #Remove first block which is corrupted
		    t,v,block_number = RemoveFirstBlock(t,v,block_number)


		    t,v,block_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,kTime,kVolt)
		    #print(t[0])
		    #print(t[0])
		    if(int(channel)==24):
		    	t=t-0.3125


		    t,v= RemoveBackwardsSamples(t,v)
		    f = Akima1DInterpolator(t,v)
		    v = f(t_univ)
		    #print(v)
		    all_volts[i,:]=v
		    #all_t[i,:]=t
		    all_blocks[i]=block_number

		    gr1.Delete()
		    usefulEvent.Delete()
		    #print(v)
		    #print(t_univ)
		    #plt.plot(t_univ,v)
		    #plt.show()

		#all_t = all_t[~np.all(all_volts==0,axis=1)]
		all_volts = all_volts[~np.all(all_volts == 0, axis=1)]
		all_blocks = all_blocks[all_blocks !=701]
		t_univ = t_univ-21.0

		np.save(global_constants.DIR_LOC+"data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",t_univ)
		np.save(global_constants.DIR_LOC+"data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_volts)
		np.save(global_constants.DIR_LOC+"data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kCalPulser)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_blocks)


		return(t_univ,all_volts)