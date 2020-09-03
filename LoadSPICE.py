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

#import BeaconTau as bt
#import beacon_data_reader
import requests
from bs4 import BeautifulSoup
from scipy.fftpack import fft
import datetime as dt
from scipy.spatial.distance import pdist
import time
#import nuphase_data_reader
#sys.path.append('../')
#import global_constants
#from model.IceModel import manytraces, LoadModel
#from load.AntennaLocations import LoadPALocations
#from load.LoadData import LoadCalPulser, LoadThermalNoise, LoadDeepPulser
#from load import nuphase_data_reader
#from process.MatchAcrossSystems import get_url_paths

from load.Calibrate import Calibrator,RemoveFirstBlock,RemoveBackwardsSamples,LoadNewPed

def LoadSPICElist():

    for i in range(5124,totalEvents):#5124
        #print(i)
        test2 = eventTree.GetEntry(i)


        usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
        gr1 = usefulEvent.getGraphFromElecChan(1)
        t_buff = gr1.GetX()
        v_buff = gr1.GetY()
        n = gr1.GetN()
        t_buff.SetSize(n)
        v_buff.SetSize(n)
        v1 = np.array(v_buff,copy=True)
        t1 = np.array(t_buff,copy=True)

        gr1 = usefulEvent.getGraphFromElecChan(9)
        t_buff = gr1.GetX()
        v_buff = gr1.GetY()
        n = gr1.GetN()
        t_buff.SetSize(n)
        v_buff.SetSize(n)
        v2 = np.array(v_buff,copy=True)
        t2 = np.array(t_buff,copy=True)

        if(np.max(v1)-np.min(v1)>500 and np.max(v2)-np.min(v2)>500):
            event_list.append(i)
            print(i)
    np.save('SPICE_eventlist.npy',np.asarray(event_list))

def LoadSpiceData(channel,kPed,kTime,kVolt):
    station = '5'
    run = '3989'
    year='2018'
    date='1225'
    length = 2048

    t_file = global_constants.DIR_LOC+"data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    v_file = global_constants.DIR_LOC+"data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    #b_file = "/home/kahughes/PA_Analysis/data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"

    if(os.path.isfile(t_file) and os.path.isfile(v_file)):
        t=np.load(t_file)
        v=np.load(v_file)
        #b = np.load(b_file)
        return(t,v)

    else:

        #ped_values = DownloadPedFile(year,run,'')
        ped_values, ped_dir = LoadNewPed(run)
        ROOT.gSystem.Load("$ARA_UTIL_INSTALL_DIR/lib/libAraEvent.so")

        
        test = ROOT.TFile.Open("/project2/avieregg/ARA_cal_data/data/root/run00"+str(run)+"/event00"+str(run)+".root")
        calibrator = ROOT.AraEventCalibrator.Instance()
        eventTree = test.Get("eventTree")
        """
        except ReferenceError:
            test = ROOT.TWebFile.Open("http://icecube:skua@convey.icecube.wisc.edu/data/exp/ARA/"+year+"/filtered/L0/ARA05/"+date+"/run"+run+"/event"+run+".root")
            calibrator = ROOT.AraEventCalibrator.Instance()
            eventTree = test.Get("eventTree")
        """
        rawEvent = ROOT.RawAtriStationEvent()

        #unixTime = HK_data.unixTime
        #print('unix time is', unixTime)
        eventTree.SetBranchAddress("event",ROOT.AddressOf(rawEvent))
        totalEvents = eventTree.GetEntries()
        print('total events:', totalEvents)
        #length = 1792

        all_volts = np.zeros([totalEvents,length])
        all_t=np.zeros([totalEvents,length])
        all_blocks=np.zeros([totalEvents])+701

        event_list = []

        event_list = np.load(global_constants.DIR_LOC+'data/SPICE_eventlist.npy')
        event_times = []
        print(channel)

        for i in event_list:#totalEvents):
            if(i%10==0):
                print(i/len(event_list))
            test2 = eventTree.GetEntry(i)
            usefulEvent = ROOT.UsefulAtriStationEvent(rawEvent,ROOT.AraCalType.kNoCalib)
            gr1 = usefulEvent.getGraphFromElecChan(channel)
            HK_data = rawEvent.unixTime
            event_times.append(HK_data)
            #t_buff = gr1.GetX()
            #v_buff = gr1.GetY()
            #n = gr1.GetN()
            #t_buff.SetSize(n)
            #v_buff.SetSize(n)
            #v = np.array(v_buff,copy=True)
            #t = np.array(t_buff,copy=True)


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
            #print(kTime,kPed,kVolt)
            t,v,block_number = Calibrator(station,t,v,block_number,str(channel),length,ped_values,kPed,kTime,kVolt)
            #print(len(v))
            if(len(v)<length):
                #print('here')
                v=np.append(v,np.zeros(length-len(v)))
            if(len(t)<length):
                t=np.append(t,np.zeros(length-len(t)))
            all_volts[i,:]=v

            all_t[i,:]=t
            all_blocks[i]=block_number

            gr1.Delete()
            usefulEvent.Delete()
            #print(t)

        all_t = all_t[~np.all(all_volts==0,axis=1)]
        all_volts = all_volts[~np.all(all_volts == 0, axis=1)]
        all_blocks = all_blocks[all_blocks !=701]

        if(kTime==1):
            all_t,all_volts= RemoveBackwardsSamples(all_t,all_volts)

        #print(t,all_volts)
        np.save(global_constants.DIR_LOC+"data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_t)
        np.save(global_constants.DIR_LOC+"data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_volts)
        np.save(global_constants.DIR_LOC+"data/SavedCalibData/blocks_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy",all_blocks)
        np.save(global_constants.DIR_LOC+"data/SPICE_event_times.npy",event_times)
        if(wantBlocks==1):
            return(all_t,all_volts,all_blocks)
        else:
            return(all_t,all_volts)