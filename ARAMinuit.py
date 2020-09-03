'''
This script use time delays measured by frequency_domain_time_delays to estimate the phase
center locations of the antennas.
'''
#import ROOT
import ctypes
import numpy as np
import itertools
import os
import sys
import csv
from iminuit import Minuit
import matplotlib.pyplot as plt
import warnings
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d, Akima1DInterpolator
import matplotlib
warnings.simplefilter(action='ignore', category=FutureWarning)
#plt.ion()
#ROOT.gSystem.Load("/home/kahughes/AraSim/libAra.so")
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

def CalculateCalTimes():
    rvals = np.linspace(20.0,100.0,160)
    zs_vals = np.linspace(-175.0,-170.0,20)
    zd_vals = np.linspace(-200.0,-140.0,200)
    prop_times = np.zeros([len(rvals),len(zs_vals),len(zd_vals)])
    LoadModel(1.326,1.78,0.0202)
    for r in range(0,len(rvals)):
        print('')
        print(r)
        for zs in range(0,len(zs_vals)):
            #print(zs)
            for zd in range(0,len(zd_vals)):
                prop_times[r,zs,zd] = trace(rvals[r],zs_vals[zs],zd_vals[zd])#ch 8 is the closer one
    np.save('/home/kahughes/PA_Analysis/data/CAL_prop_times.npy',prop_times)
    return(prop_times)

def LoadCalTimes():
    prop_times = np.load('/home/kahughes/PA_Analysis/data/CAL_prop_times.npy')
    rvals = np.linspace(20.0,100.0,160)
    zs_vals = np.linspace(-175.0,-170.0,20)
    zd_vals = np.linspace(-200.0,-140.0,200)
    interpolator = RegularGridInterpolator((rvals,zs_vals,zd_vals),prop_times)
    #print(interpolator([4160,-1090.,-180]))
    return(interpolator)

def CalculatePropTimes():
    rvals = np.linspace(4100.0,4190.0,90)
    zs_vals = np.linspace(-1450,-1050,20)
    zd_vals = np.linspace(-200.0,-140.0,120)
    prop_times = np.zeros([len(rvals),len(zs_vals),len(zd_vals)])
    LoadModel(1.326,1.78,0.0202)
    for r in range(0,len(rvals)):
        #print('')
        #print(r)
        for zs in range(0,len(zs_vals)):
            #print(zs)
            for zd in range(0,len(zd_vals)):
                prop_times[r,zs,zd] = trace(rvals[r],zs_vals[zs],zd_vals[zd])#ch 8 is the closer one
    np.save('/home/kahughes/PA_Analysis/data/SPICE_prop_times.npy',prop_times)
    return(prop_times)

def LoadPropTimes():
    prop_times = np.load('/home/kahughes/PA_Analysis/data/SPICE_prop_times.npy')
    rvals = np.linspace(4100.0,4190.0,90)
    zs_vals = np.linspace(-1450.0,-1050.0,20)
    zd_vals = np.linspace(-200.0,-140.0,120)
    interpolator = RegularGridInterpolator((rvals,zs_vals,zd_vals),prop_times)
    #print(interpolator([4160,-1090.,-180]))
    return(interpolator)

def LoadLocalSPICE(run,channel,kPed,kTime,kVolt):
    t_file = "/home/kahughes/PA_Analysis/data/SavedCalibData/time_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    v_file = "/home/kahughes/PA_Analysis/data/SavedCalibData/volts_"+str(run)+"_"+str(channel)+str(kPed)+str(kTime)+str(kVolt)+".npy"
    t=np.load(t_file)
    v=np.load(v_file)
    return(t,v)
def oldvals():
    """
    ant1_pxb = +43.225626
    ant1_pyb = +51.534588
    ant1_pzb = -192.607929
    ant1_pxt = +43.209428
    ant1_pyt = +51.488248
    ant1_pzt = -162.990586

    ant2_pxb = +21.430105
    ant2_pyb = +63.319284
    ant2_pzb = -194.037016
    ant2_pxt = +21.405678
    ant2_pyt = +63.278645
    ant2_pzt = -164.393729

    ant3_pxb = +27.760625
    ant3_pyb = +22.522423
    ant3_pzb = -189.013727
    ant3_pxt = +27.856493
    ant3_pyt = +22.617241
    ant3_pzt = -159.217541

    ant4_pxb = +5.391217
    ant4_pyb = +39.633913
    ant4_pzb = -175.450046
    ant4_pxt = +5.389867
    ant4_pyt = +39.669012
    ant4_pzt = -144.942933


    cal_pz = -171.278558
    SPICE_x_guess = +3072.295822
    SPICE_y_guess = +2861.7688200

    c0p = +138.485012
    c1p = +17.547130
    c8p = +133.911865
    c9p = +18.789163
    c16p = +131.291872
    c17p = +13.460643
    c24p = +132.996322
    c25p = +16.938025

    #ORIGINAL GUESSES
    ant2_pxb = 19.79
    ant2_pyb = 64.73
    ant2_pzb = -195.31
    ant2_pxt = 19.79
    ant2_pyt = 64.73
    ant2_pzt = -165.46

    ant1_pxb = 44.89
    ant1_pyb = 49.59
    ant1_pzb = -194.52
    ant1_pxt = 44.89
    ant1_pyt = 49.59
    ant1_pzt = -164.67

    ant3_pxb = 27.06
    ant3_pyb = 21.85
    ant3_pzb = -191.28
    ant3_pxt = 27.06
    ant3_pyt = 21.85
    ant3_pzt = -161.34

    ant4_pxb = 5.27
    ant4_pyb = 38.95
    ant4_pzb = -177.72
    ant4_pxt = 5.27
    ant4_pyt = 38.95
    ant4_pzt = -145.01

    cal_pz = -173.28
    SPICE_x_guess = +3072.295822
    SPICE_y_guess = +2861.7688200
    c0p = 132.3
    c1p = 18.62
    c8p = 132.3
    c9p = 18.62#32.93#
    c16p = 132.3#125.85#
    c17p = 18.62#27.79#
    c24p = 132.3#129.07#
    c25p = 18.62#12.51#

    #FIRST FIT GUESSES WITH BAD PEDESTALS
    ant1_pxb= 44.75447143486855
    ant1_pyb= 49.534417695859176
    ant1_pzb= -193.3476220529692
    ant1_pxt= 44.75447143486855
    ant1_pyt= 49.534417695859176
    ant1_pzt= -163.57191431085786

    ant2_pxb= 19.627862410578622
    ant2_pyb= 63.99912124281646
    ant2_pzb= -194.7726712019533
    ant2_pxt= 19.627862410578622
    ant2_pyt= 63.99912124281646
    ant2_pzt= -164.82568829824638

    ant3_pxb= 27.602174341653324
    ant3_pyb= 22.383567222629054
    ant3_pzb= -189.29760854752965
    ant3_pxt= 27.602174341653324
    ant3_pyt= 22.383567222629054
    ant3_pzt= -159.73247856193603

    ant4_pxb= 5.60974680169312
    ant4_pyb= 39.091309722809065
    ant4_pzb= -176.15729903644115
    ant4_pxt= 5.60974680169312
    ant4_pyt= 39.091309722809065
    ant4_pzt= -145.56616424530017
    cal_pz= -171.92011053883772
    SPICE_x_guess= 3068.7831991294565
    SPICE_y_guess= 2866.73590803123
    c0p= 139.484743343324
    c1p= 18.554970476620106
    c8p= 131.52293639009153
    c9p= 16.2403487932527
    c16p=132.573129390966
    c17p= 14.436596497767997
    c24p= 133.4764217409696
    c25p= 16.93045009201646
    

    ant1_pxb= 44.75941265639587
    ant1_pyb= 49.53904655592454
    ant1_pzb= -193.4326087997942
    ant1_pxt= 44.75941265639587
    ant1_pyt= 49.53904655592454
    ant1_pzt= -163.67166420392869
    ant2_pxb= 19.63557266551782
    ant2_pyb= 64.00678139537834
    ant2_pzb= -194.84463411215245
    ant2_pxt= 19.63557266551782
    ant2_pyt= 64.00678139537834
    ant2_pzt= -164.87955697468567
    ant3_pxb= 27.604821294543445
    ant3_pyb= 22.38542592546542
    ant3_pzb= -189.32053259158548
    ant3_pxt= 27.604821294543445
    ant3_pyt= 22.38542592546542
    ant3_pzt= -159.7995021052749
    ant4_pxb= 5.606121325342855
    ant4_pyb= 39.111087498262094
    ant4_pzb= -176.1951057304098
    ant4_pxt= 5.606121325342855
    ant4_pyt= 39.11108749826209
    ant4_pzt= -145.69465701904616
    cal_pz= -172.0247301502089
    SPICE_x_guess= 3067.419184577016
    SPICE_y_guess= 2865.050201559079
    c0p= 139.48859918384522
    c1p= 18.520780204100898
    c8p= 131.53828670226974
    c9p= 16.16398424839889
    c16p= 132.51554135576768
    c17p= 14.429473595270691
    c24p= 133.6856809759707
    c25p= 16.84802953929185
    """
    #plt.figure(0)
    #plt.scatter([ant1_px,ant2_px,ant3_px,ant4_px],[ant1_py,ant2_py,ant3_py,ant4_py])
    #plt.scatter(0.0,0.0)
    #plt.scatter(SPICE_x_guess,SPICE_y_guess)
    #plt.show()

    #print('try sample f:')
    #init_chi2 = f(ant1_pxb,ant1_pyb,ant1_pzb,ant1_pxt,ant1_pyt,ant1_pzt, ant2_px, ant2_py, ant2_pzb, ant2_pzt,ant3_px,ant3_py,ant3_pzb,ant3_pzt,ant4_px,ant4_py,ant4_pzb,ant4_pzt,cal_pz,SPICE_x_guess,SPICE_y_guess,c0p,c1p,c8p,c9p,c16p,c17p,c24p,c25p)
    #print('Inital chi2 is:', init_chi2)
def LoadSPICEwf(pairs):

    t0,v0 = LoadLocalSPICE("3989",int(0),1,1,1)
    t1,v1 = LoadLocalSPICE("3989",int(1),1,1,1)
    t8,v8 = LoadLocalSPICE("3989",int(8),1,1,1)
    t9,v9 = LoadLocalSPICE("3989",int(9),1,1,1)
    t16,v16 = LoadLocalSPICE("3989",int(16),1,1,1)
    t17,v17 = LoadLocalSPICE("3989",int(17),1,1,1)
    t24,v24 = LoadLocalSPICE("3989",int(24),1,0,0)
    #t24 = t24-t24[0,0]
    t24 = t24-(t24[0,0]-20.0)
    t25,v25 = LoadLocalSPICE("3989",int(25),1,1,1)
    delay_dict_np = np.load('delay_dict.npy',allow_pickle=True,encoding='latin1')
    depth_dict_np = np.load('depth_dict.npy',allow_pickle=True,encoding='latin1')
    delays = {}
    depths = {}
    for i in pairs:
        delays.update({i:delay_dict_np.item().get(i)})
        depths.update({i:depth_dict_np.item().get(i)})
    #print('here are all depths', depths[(0,16)])
    SPICE_depths = np.linspace(-1050,-1450,20)
    interp_delays = {}
    for i in pairs:
        #i=(16,17)
        #if(i[0]==0 or i[0]==1 or i[0]==8 or i[0]==9):
        #    continue
        print(i)
        top_vs = [0,8,24]
        bottom_vs = [1,9,25]
        if(i[0]==1 and i[1]==9):
            cut_vals = 2.0
            linear_fit = 1
        elif(i[0]==1 and i[1]==17):
            cut_vals = 5.0
            linear_fit = 0
        elif(i[0]==9 and i[1]==17):
            cut_vals = 9.0
            linear_fit = 0
        elif(i[0]==25 and i[1]==17):
            cut_vals = 7.0
            linear_fit = 0
        elif(i[0]==0 and i[1]==25):
            cut_vals=10.0
            linear_fit = 0
        elif(i[0]==0 and i[1]==16):
            cut_vals=3.0
            linear_fit = 0

        elif((i[0] in top_vs and i[1] in bottom_vs) or(i[0] in bottom_vs and i[1] in top_vs)):
            cut_vals = 15.0
            linear_fit = 0
        elif(i[0]==16 or i[1]==16):
            cut_vals = 15.0
            linear_fit = 0
        elif(i[0]==17 or i[1]==17):
            cut_vals = 25.0
            linear_fit = 0

        #elif(i[0]==0 and i[1]==24):
        #    cut_vals = 3.0
        #
        else:
            cut_vals = 2.0
            linear_fit = 1

        this_delay = np.asarray(delays[i])
        this_depth = np.asarray(depths[i])
        print(this_delay,this_depth,cut_vals,linear_fit)
        points_t = []
        if(i[0]==1 and i[1]==17):
            this_mean=-222.5
            points_t = [-218.6,-226.506]
            points_d = [-990,-1447.72]
        #elif(i[0]==8 and i[1]==17):
        #    this_mean=-105.0
        #elif(i[0]==9 and i[1]==17):
        #    this_mean=-178.0
        elif(i[0]==25 and i[1]==17):
            this_mean=-42.0
        elif(i[0]==0 and i[1]==8):
            this_mean=-40.0
        elif(i[0]==0 and i[1]==24):
            this_mean=-178.0
        elif(i[0]==0 and i[1]==25):
            this_mean=-100.0
            points_t = [-90.4,-106.82]
            points_d = [-1029,-1450]
        elif(i[0]==0 and i[1]==16):
            this_mean = np.mean(this_delay)
            points_t = [-219.27,-222.7]
            points_d = [-1075,-1336]
        elif(i[0]==0 and i[1]==17):
            this_mean = np.mean(this_delay)
            points_t = [-153,-124.485]
            points_d = [-1445,-997.0]
        elif(i[0]==1 and i[1]==9):
            this_mean = -45.0
        elif(i[0]==1 and i[1]==24):
            this_mean = -260.0
            points_t = [-272.5,-250]
            points_d = [-985,-1450]
        elif(i[0]==1 and i[1]==25):
            this_mean=-181.5
        elif(i[0]==1 and i[1]==16):
            this_mean=-303.0
            points_t = [-311.0,-297.3]
            points_d = [-995.0,-1450.0]
        elif(i[0]==8 and i[1]==25):
            this_mean=-60.0
            points_t=[-48.4,-68.53]
            points_d=[-1033.8,-1454.8]
        elif(i[0]==8 and i[1]==16):
            this_mean = np.mean(this_delay)
            points_t=[-177.25,-185.309]
            points_d=[-988.6,-1456.26]
        elif(i[0]==8 and i[1]==17):
            this_mean = np.mean(this_delay)
            points_t=[-83.9,-113.912]
            points_d=[-984.23,-1454.34]
        elif(i[0]==9 and i[1]==25):
            this_mean=-135.0
        elif(i[0]==9 and i[1]==16):
            this_mean = np.mean(this_delay)
            points_t=[-264.225,-252.25]
            points_d=[-990.0,-1461.0]
        elif(i[0]==9 and i[1]==17):
            this_mean = -178.0
            points_t=[-171.37,-181.68]
            points_d=[-990.2,-1450.5]
        elif(i[0]==16 and i[1]==17):
            this_mean=80.0
            points_t=[93.42,70.63]
            points_d=[-987.1,-1460.26]
        else:
            this_mean = np.mean(this_delay)

        print(this_mean)
        #find close values to mean:
        my_vals = np.where(np.abs(this_delay-this_mean)<cut_vals)
        #print(my_vals)
        #plt.scatter(this_delay,this_depth)
        #plt.plot(new_delay,new_depth,color='darkblue',lw=2.5)
        #plt.show()
        this_depth = this_depth[my_vals]
        this_delay = this_delay[my_vals]
        #plt.scatter(this_delay,this_depth)
        #plt.show()
        if(len(points_t)>0):
            fitvals = np.polyfit(points_t,points_d,1)
        else:

            fitvals = np.polyfit(this_delay,this_depth,1)
        test = np.abs(this_depth-(fitvals[1]+fitvals[0]*this_delay))
        print(test)
        """
        plt.figure(0,facecolor='w')

        plt.scatter(this_delay[np.where(this_depth>-1470)],this_depth[np.where(this_depth>-1470)],label='All Data')
        """
        #plt.plot(this_delay,(fitvals[1]+fitvals[0]*this_delay),color='red')
        #plt.plot(new_delay,new_depth,color='darkblue',lw=2.5)
        #plt.show()
        if(linear_fit==0):
            if((i[0]==9 or i[0]==16) and i[1]==17):
                this_delay = this_delay[test<50.0]
                this_depth = this_depth[test<50.0]
            else:
                this_delay = this_delay[test<50.0]
                this_depth = this_depth[test<50.0]

        #print('here')
        new_delay = []
        new_depth = []
        for j in range(0,int(len(this_delay)/10)):
            new_delay.append(np.mean(this_delay[j*10:j*10+10]))
            new_depth.append(np.mean(this_depth[j*10:j*10+10]))
        new_delay = np.asarray(new_delay)
        new_depth = np.asarray(new_depth)
        #print('here')
        #plt.figure(0,facecolor='w')
        #plt.scatter(this_delay,this_depth)
        #plt.plot(new_delay,new_depth,color='darkblue',lw=2.5)
        #plt.show()

        f_int = interp1d(new_depth,new_delay)
        print(np.min(new_depth),np.max(new_depth))
        print(SPICE_depths)
        temp_delay = np.zeros(len(SPICE_depths))
        for d in range(0,len(SPICE_depths)):
            temp_delay[d]=f_int(SPICE_depths[d])
        print(SPICE_depths)
        """
        plt.scatter(temp_delay,SPICE_depths,color='red',label='Points Used for Optimization')
        plt.grid()
        plt.legend()
        plt.xlabel('Raw Time Delay (ns)')
        plt.ylabel('SPIceCore Depth (m)')
        plt.show()
        """
        #plt.figure(0,facecolor='w')
        #plt.scatter(this_delay,this_depth)
        #plt.plot(temp_delay,SPICE_depths)
        #plt.show()
        interp_delays.update({i:np.asarray(temp_delay)})
    return(interp_delays,SPICE_depths)

def FindR(d1x,d1y,d2x,d2y):
    return(np.sqrt((d1y-d2y)**2+(d1x-d2x)**2))

def f2(ant1_xb, ant1_yb, ant1_zb, ant1_xt, ant1_yt, ant1_zt, ant2_xb, ant2_yb, ant2_zb, ant2_xt, ant2_yt, ant2_zt, ant3_xb, ant3_yb, ant3_zb, ant3_xt, ant3_yt,  ant3_zt, ant4_xb, ant4_yb, ant4_zb, ant4_xt, ant4_yt,ant4_zt, cal_z, SPICE_x,SPICE_y,c0,c1,c8,c9,c16,c17,c24,c25):
    #To generalize, look into from_array_func Minuit initializer.


    try:

        #fixing the locations of antenna zero.
        coord_dict = {0:[ant1_xb,ant1_yb,ant1_zb],1:[ant1_xt,ant1_yt,ant1_zt],8:[ant2_xb,ant2_yb,ant2_zb],9:[ant2_xt,ant2_yt,ant2_zt],16:[ant4_xb,ant4_yb,ant4_zb],17:[ant4_xt,ant4_yt,ant4_zt],24:[ant3_xb,ant3_yb,ant3_zb],25:[ant3_xt,ant3_yt,ant3_zt]}
        #coord_dict = {0:[ant1_xb,ant1_yb,ant1_zb],1:[ant1_xb,ant1_yb,ant1_zt],8:[ant2_xb,ant2_yb,ant2_zb],9:[ant2_xb,ant2_yb,ant2_zt],16:[ant4_xb,ant4_yb,ant4_zb],17:[ant4_xb,ant4_yb,ant4_zt],24:[ant3_xb,ant3_yb,ant3_zb],25:[ant3_xb,ant3_yb,ant3_zt]}

        cable_delays = {0:c0,1:c1,8:c8,9:c9,16:c16,17:c17,24:c24,25:c25}

        SPICE_r_dict = {}
        cal_r_dict = {}
        for j in [0,1,8,9,16,17,24,25]:
            #print(j)
            #print(FindR(coord_dict[j][0],coord_dict[j][1],SPICE_x,SPICE_y))
            SPICE_r_dict.update({j:FindR(coord_dict[j][0],coord_dict[j][1],SPICE_x,SPICE_y)})
            cal_r_dict.update({j:FindR(coord_dict[j][0],coord_dict[j][1],0.0,0.0)})
        #print(SPICE_r_dict)
        #print(ant1_zb,ant1_zt,ant2_zb,ant2_zt,ant3_zb,ant3_zt,ant4_zb,ant4_zt)
        #print(cal_r_dict)
        chi_2 = 0.0
        for i in pairs:
            counter = 0
            #print(i[0],SPICE_r_dict[i[0]],coord_dict[i[0]][2])

            #if(i[1]==8):
            #    print('')
            #    print(ant2_x,ant2_y)
            #    print(coord_dict[8][0],coord_dict[8][1])
            #    print(i[1],SPICE_r_dict[i[1]],coord_dict[i[1]][2])
            #print(cal_r_dict[i[0]],cal_r_dict[i[1]],coord_dict[i[0]][2],coord_dict[i[1]][2])
            t1 = cal_times([cal_r_dict[i[0]],cal_z,coord_dict[i[0]][2]])
            t2 = cal_times([cal_r_dict[i[1]],cal_z,coord_dict[i[1]][2]])
            #print('')
            #print((t2-t1)[0],cable_delays[i[0]],cable_delays[i[1]],cal_delays[i])
            chi_2 = chi_2+(t2-t1-cable_delays[i[0]]+cable_delays[i[1]]+cal_delays[i])**2

            for zs in SPICE_depths:
                #print((t2-t1)[0],cable_delays[i[0]],cable_delays[i[1]],int_delays[i][counter])

                t1 = interp([SPICE_r_dict[i[0]],zs,coord_dict[i[0]][2]])
                t2 = interp([SPICE_r_dict[i[1]],zs,coord_dict[i[1]][2]])
                #print((t2-t1)[0],cable_delays[i[0]],cable_delays[i[1]],int_delays[i][counter])
                chi_2 = chi_2+(t2-t1-cable_delays[i[0]]+cable_delays[i[1]]+int_delays[i][counter])**2
                counter = counter + 1

        print(chi_2)
        return chi_2
    except Exception as e:
        print('Error in f')
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

def f1(ant1_xb, ant1_yb, ant1_zb, ant1_zt, ant2_xb, ant2_yb, ant2_zb, ant2_zt, ant3_xb, ant3_yb, ant3_zb,  ant3_zt, ant4_xb, ant4_yb, ant4_zb, ant4_zt, cal_z, SPICE_x,SPICE_y,c0,c1,c8,c9,c16,c17,c24,c25):
    #To generalize, look into from_array_func Minuit initializer.

    try:

        #fixing the locations of antenna zero.
        #coord_dict = {0:[ant1_xb,ant1_yb,ant1_zb],1:[ant1_xt,ant1_yt,ant1_zt],8:[ant2_xb,ant2_yb,ant2_zb],9:[ant2_xt,ant2_yt,ant2_zt],16:[ant4_xb,ant4_yb,ant4_zb],17:[ant4_xt,ant4_yt,ant4_zt],24:[ant3_xb,ant3_yb,ant3_zb],25:[ant3_xt,ant3_yt,ant3_zt]}
        coord_dict = {0:[ant1_xb,ant1_yb,ant1_zb],1:[ant1_xb,ant1_yb,ant1_zt],8:[ant2_xb,ant2_yb,ant2_zb],9:[ant2_xb,ant2_yb,ant2_zt],16:[ant4_xb,ant4_yb,ant4_zb],17:[ant4_xb,ant4_yb,ant4_zt],24:[ant3_xb,ant3_yb,ant3_zb],25:[ant3_xb,ant3_yb,ant3_zt]}

        cable_delays = {0:c0,1:c1,8:c8,9:c9,16:c16,17:c17,24:c24,25:c25}

        SPICE_r_dict = {}
        cal_r_dict = {}
        for j in [0,1,8,9,16,17,24,25]:
            #print(j)
            #print(FindR(coord_dict[j][0],coord_dict[j][1],SPICE_x,SPICE_y))
            SPICE_r_dict.update({j:FindR(coord_dict[j][0],coord_dict[j][1],SPICE_x,SPICE_y)})
            cal_r_dict.update({j:FindR(coord_dict[j][0],coord_dict[j][1],0.0,0.0)})
        #print(SPICE_r_dict)
        #print(ant1_zb,ant1_zt,ant2_zb,ant2_zt,ant3_zb,ant3_zt,ant4_zb,ant4_zt)
        #print(cal_r_dict)
        chi_2 = 0.0
        for i in pairs:
            counter = 0
            #print(i[0],SPICE_r_dict[i[0]],coord_dict[i[0]][2])

            #if(i[1]==8):
            #    print('')
            #    print(ant2_x,ant2_y)
            #    print(coord_dict[8][0],coord_dict[8][1])
            #    print(i[1],SPICE_r_dict[i[1]],coord_dict[i[1]][2])
            #print(cal_r_dict[i[0]],cal_r_dict[i[1]],coord_dict[i[0]][2],coord_dict[i[1]][2])
            t1 = cal_times([cal_r_dict[i[0]],cal_z,coord_dict[i[0]][2]])
            t2 = cal_times([cal_r_dict[i[1]],cal_z,coord_dict[i[1]][2]])
            #print('')
            #print((t2-t1)[0],cable_delays[i[0]],cable_delays[i[1]],cal_delays[i])
            chi_2 = chi_2+(t2-t1-cable_delays[i[0]]+cable_delays[i[1]]+cal_delays[i])**2

            for zs in SPICE_depths:
                #print((t2-t1)[0],cable_delays[i[0]],cable_delays[i[1]],int_delays[i][counter])

                t1 = interp([SPICE_r_dict[i[0]],zs,coord_dict[i[0]][2]])
                t2 = interp([SPICE_r_dict[i[1]],zs,coord_dict[i[1]][2]])
                #print((t2-t1)[0],cable_delays[i[0]],cable_delays[i[1]],int_delays[i][counter])
                chi_2 = chi_2+(t2-t1-cable_delays[i[0]]+cable_delays[i[1]]+int_delays[i][counter])**2
                counter = counter + 1

        print(chi_2)
        return chi_2
    except Exception as e:
        print('Error in f')
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

def FirstIteration():
    print('done')
    
    #SPICE_delays = np.load('all_SPICE_delays_int.npy',allow_pickle=True)
    #SPICE_depths = np.load('all_SPICE_depths_int.npy',allow_pickle=True)
    #print('SPICE is', SPICE_delays,SPICE_depths)

    ant2_pxb = 19.79
    ant2_pyb = 64.73
    ant2_pzb = -195.31
    ant2_pxt = 19.79
    ant2_pyt = 64.73
    ant2_pzt = -165.46

    ant1_pxb = 44.89
    ant1_pyb = 49.59
    ant1_pzb = -194.52
    ant1_pxt = 44.89
    ant1_pyt = 49.59
    ant1_pzt = -164.67

    ant3_pxb = 27.06
    ant3_pyb = 21.85
    ant3_pzb = -191.28
    ant3_pxt = 27.06
    ant3_pyt = 21.85
    ant3_pzt = -161.34

    ant4_pxb = 5.27
    ant4_pyb = 38.95
    ant4_pzb = -177.72
    ant4_pxt = 5.27
    ant4_pyt = 38.95
    ant4_pzt = -145.01


    cal_pz = -171.278558
    SPICE_x_guess = +3072.295822
    SPICE_y_guess = +2861.7688200

    c0p = +138.485012
    c1p = +17.547130
    c8p = +133.911865
    c9p = +18.789163
    c16p = +131.291872
    c17p = +13.460643
    c24p = +132.996322
    c25p = +16.938025



    rvals = np.linspace(4100.0,4190.0,90)
    zs_vals = np.linspace(-1450,-1050,20)
    zd_vals = np.linspace(-200.0,-140.0,120)

    initial_step = 0.1 #m
    print('trying minimization now')
    #-12 ft on pulser locations relative to antennas to account for additional mast elevation.
    m = Minuit(     f1,\
                    ant1_xb=ant1_pxb,\
                    ant1_yb=ant1_pyb,\
                    ant1_zb=ant1_pzb,\
                    ant1_zt=ant1_pzt,\
                    ant2_xb=ant2_pxb,\
                    ant2_yb=ant2_pyb,\
                    ant2_zb=ant2_pzb,\
                    ant2_zt=ant2_pzt,\
                    ant3_xb=ant3_pxb,\
                    ant3_yb=ant3_pyb,\
                    ant3_zb=ant3_pzb,\
                    ant3_zt=ant3_pzt,\
                    ant4_xb=ant4_pxb,\
                    ant4_yb=ant4_pyb,\
                    ant4_zb=ant4_pzb,\
                    ant4_zt=ant4_pzt,\
                    cal_z = cal_pz,\
                    SPICE_x=SPICE_x_guess,\
                    SPICE_y=SPICE_y_guess,\
                    c0 = c0p,\
                    c1 = c1p,\
                    c8 = c8p,\
                    c9 = c9p,\
                    c16 = c16p,\
                    c17 = c17p,\
                    c24 = c24p,\
                    c25 = c25p,\
                    error_ant1_xb=initial_step,\
                    error_ant1_yb=initial_step,\
                    error_ant1_zb=initial_step,\
                    error_ant1_zt=initial_step,\
                    error_ant2_xb=initial_step,\
                    error_ant2_yb=initial_step,\
                    error_ant2_zb=initial_step,\
                    error_ant2_zt=initial_step,\
                    error_ant3_xb=initial_step,\
                    error_ant3_yb=initial_step,\
                    error_ant3_zb=initial_step,\
                    error_ant3_zt=initial_step,\
                    error_ant4_xb=initial_step,\
                    error_ant4_yb=initial_step,\
                    error_ant4_zb=initial_step,\
                    error_ant4_zt=initial_step,\
                    error_cal_z=initial_step,\
                    error_SPICE_x=initial_step,\
                    error_SPICE_y=initial_step,\
                    error_c0 = initial_step,\
                    error_c1 = initial_step,\
                    error_c8 = initial_step,\
                    error_c9 = initial_step,\
                    error_c16 = initial_step,\
                    error_c17 = initial_step,\
                    error_c24 = initial_step,\
                    error_c25 = initial_step,\
                    errordef = 1.0,\
                    limit_ant1_xb=(ant1_pxb-10,ant1_pxb+10),\
                    limit_ant1_yb=(ant1_pyb-10,ant1_pyb+10),\
                    limit_ant1_zb=(ant1_pzb-1,ant1_pzb+2),\
                    limit_ant1_zt=(ant1_pzt-1,ant1_pzt+2),\
                    limit_ant2_xb=(ant2_pxb-10,ant2_pxb+10),\
                    limit_ant2_yb=(ant2_pyb-10,ant2_pyb+10),\
                    limit_ant2_zb=(ant2_pzb-1,ant2_pzb+2),\
                    limit_ant2_zt=(ant2_pzt-1,ant2_pzt+2),\
                    limit_ant3_xb=(ant3_pxb-10,ant3_pxb+10),\
                    limit_ant3_yb=(ant3_pyb-10,ant3_pyb+10),\
                    limit_ant3_zb=(ant3_pzb-1,ant3_pzb+2),\
                    limit_ant3_zt=(ant3_pzt-1,ant3_pzt+2),\
                    limit_ant4_xb=(ant4_pxb-10,ant4_pxb+10),\
                    limit_ant4_yb=(ant4_pyb-10,ant4_pyb+10),\
                    limit_ant4_zb=(ant4_pzb-1,ant4_pzb+2),\
                    limit_ant4_zt=(ant4_pzt-1,ant4_pzt+2),\
                    limit_SPICE_x=(SPICE_x_guess-5.0,SPICE_x_guess+5.0),\
                    limit_SPICE_y=(SPICE_y_guess-5.0,SPICE_y_guess+5.0),\
                    limit_c0 = (c0p-15.0,c0p+15.0),\
                    limit_c1 = (c1p-15.0,c1p+15.0),\
                    limit_c8 = (c8p-15.0,c8p+15.0),\
                    limit_c9 = (c9p-15.0,c9p+15.0),\
                    limit_c16 = (c16p-15.0,c16p+15.0),\
                    limit_c17 = (c17p-15.0,c17p+15.0),\
                    limit_c24 = (c24p-15.0,c24p+15.0),\
                    limit_c25 = (c25p-15.0,c25p+15.0),\
                    limit_cal_z = (-175.0,-170.0)
                    )

    result = m.migrad()
    values = m.values
    print(result)
    print(m.values)
    #m.hesse()
    #m.minos()
    #print(m.get_fmin())
    print(result)

    #print(m.values['ant1_xb'])
    print(m.values[0])
    np.save('first_iteration.npy',m.values)
    return(values)

def SecondIteration(values):
    ant1_pxb = values[0]
    ant1_pyb = values[1]
    ant1_pzb = values[2]
    ant1_pxt = values[0]
    ant1_pyt = values[1]
    ant1_pzt = values[3]

    ant2_pxb = values[4]
    ant2_pyb = values[5]
    ant2_pzb = values[6]
    ant2_pxt = values[4]
    ant2_pyt = values[5]
    ant2_pzt = values[7]

    ant3_pxb = values[8]
    ant3_pyb = values[9]
    ant3_pzb = values[10]
    ant3_pxt = values[8]
    ant3_pyt = values[9]
    ant3_pzt = values[11]

    ant4_pxb = values[12]
    ant4_pyb = values[13]
    ant4_pzb = values[14]
    ant4_pxt = values[12]
    ant4_pyt = values[13]
    ant4_pzt = values[15]


    cal_pz = values[16]
    SPICE_x_guess = values[17]
    SPICE_y_guess = values[18]

    c0p = values[19]
    c1p = values[20]
    c8p = values[21]
    c9p = values[22]
    c16p = values[23]
    c17p = values[24]
    c24p = values[25]
    c25p = values[26]
    initial_step = 0.1 #m
    m = Minuit(     f2,\
                    ant1_xb=ant1_pxb,\
                    ant1_yb=ant1_pyb,\
                    ant1_zb=ant1_pzb,\
                    ant1_xt=ant1_pxt,\
                    ant1_yt=ant1_pyt,\
                    ant1_zt=ant1_pzt,\
                    ant2_xb=ant2_pxb,\
                    ant2_yb=ant2_pyb,\
                    ant2_zb=ant2_pzb,\
                    ant2_xt=ant2_pxt,\
                    ant2_yt=ant2_pyt,\
                    ant2_zt=ant2_pzt,\
                    ant3_xb=ant3_pxb,\
                    ant3_yb=ant3_pyb,\
                    ant3_zb=ant3_pzb,\
                    ant3_xt=ant3_pxt,\
                    ant3_yt=ant3_pyt,\
                    ant3_zt=ant3_pzt,\
                    ant4_xb=ant4_pxb,\
                    ant4_yb=ant4_pyb,\
                    ant4_zb=ant4_pzb,\
                    ant4_xt=ant4_pxt,\
                    ant4_yt=ant4_pyt,\
                    ant4_zt=ant4_pzt,\
                    cal_z = cal_pz,\
                    SPICE_x=SPICE_x_guess,\
                    SPICE_y=SPICE_y_guess,\
                    c0 = c0p,\
                    c1 = c1p,\
                    c8 = c8p,\
                    c9 = c9p,\
                    c16 = c16p,\
                    c17 = c17p,\
                    c24 = c24p,\
                    c25 = c25p,\
                    error_ant1_xb=initial_step,\
                    error_ant1_yb=initial_step,\
                    error_ant1_zb=initial_step,\
                    error_ant1_xt=initial_step,\
                    error_ant1_yt=initial_step,\
                    error_ant1_zt=initial_step,\
                    error_ant2_xb=initial_step,\
                    error_ant2_yb=initial_step,\
                    error_ant2_zb=initial_step,\
                    error_ant2_xt=initial_step,\
                    error_ant2_yt=initial_step,\
                    error_ant2_zt=initial_step,\
                    error_ant3_xb=initial_step,\
                    error_ant3_yb=initial_step,\
                    error_ant3_zb=initial_step,\
                    error_ant3_xt=initial_step,\
                    error_ant3_yt=initial_step,\
                    error_ant3_zt=initial_step,\
                    error_ant4_xb=initial_step,\
                    error_ant4_yb=initial_step,\
                    error_ant4_zb=initial_step,\
                    error_ant4_xt=initial_step,\
                    error_ant4_yt=initial_step,\
                    error_ant4_zt=initial_step,\
                    error_cal_z=initial_step,\
                    error_SPICE_x=initial_step,\
                    error_SPICE_y=initial_step,\
                    error_c0 = initial_step,\
                    error_c1 = initial_step,\
                    error_c8 = initial_step,\
                    error_c9 = initial_step,\
                    error_c16 = initial_step,\
                    error_c17 = initial_step,\
                    error_c24 = initial_step,\
                    error_c25 = initial_step,\
                    errordef = 1.0,\
                    limit_ant1_xb=(ant1_pxb-10,ant1_pxb+10),\
                    limit_ant1_yb=(ant1_pyb-10,ant1_pyb+10),\
                    limit_ant1_zb=(ant1_pzb-1,ant1_pzb+2),\
                    limit_ant1_xt=(ant1_pxt-10,ant1_pxt+10),\
                    limit_ant1_yt=(ant1_pyt-10,ant1_pyt+10),\
                    limit_ant1_zt=(ant1_pzt-1,ant1_pzt+2),\
                    limit_ant2_xb=(ant2_pxb-10,ant2_pxb+10),\
                    limit_ant2_yb=(ant2_pyb-10,ant2_pyb+10),\
                    limit_ant2_zb=(ant2_pzb-1,ant2_pzb+2),\
                    limit_ant2_xt=(ant2_pxt-10,ant2_pxt+10),\
                    limit_ant2_yt=(ant2_pyt-10,ant2_pyt+10),\
                    limit_ant2_zt=(ant2_pzt-1,ant2_pzt+2),\
                    limit_ant3_xb=(ant3_pxb-10,ant3_pxb+10),\
                    limit_ant3_yb=(ant3_pyb-10,ant3_pyb+10),\
                    limit_ant3_zb=(ant3_pzb-1,ant3_pzb+2),\
                    limit_ant3_xt=(ant3_pxt-10,ant3_pxt+10),\
                    limit_ant3_yt=(ant3_pyt-10,ant3_pyt+10),\
                    limit_ant3_zt=(ant3_pzt-1,ant3_pzt+2),\
                    limit_ant4_xb=(ant4_pxb-10,ant4_pxb+10),\
                    limit_ant4_yb=(ant4_pyb-10,ant4_pyb+10),\
                    limit_ant4_zb=(ant4_pzb-1,ant4_pzb+2),\
                    limit_ant4_xt=(ant4_pxt-10,ant4_pxt+10),\
                    limit_ant4_yt=(ant4_pyt-10,ant4_pyt+10),\
                    limit_ant4_zt=(ant4_pzt-1,ant4_pzt+2),\
                    limit_SPICE_x=(SPICE_x_guess-5.0,SPICE_x_guess+5.0),\
                    limit_SPICE_y=(SPICE_y_guess-5.0,SPICE_y_guess+5.0),\
                    limit_c0 = (c0p-15.0,c0p+15.0),\
                    limit_c1 = (c1p-15.0,c1p+15.0),\
                    limit_c8 = (c8p-15.0,c8p+15.0),\
                    limit_c9 = (c9p-15.0,c9p+15.0),\
                    limit_c16 = (c16p-15.0,c16p+15.0),\
                    limit_c17 = (c17p-15.0,c17p+15.0),\
                    limit_c24 = (c24p-15.0,c24p+15.0),\
                    limit_c25 = (c25p-15.0,c25p+15.0),\
                    limit_cal_z = (-175.0,-170.0)
                    )

    result = m.migrad()
    print(result)
    print(m.values)
    #m.hesse()
    #m.minos()
    #print(m.get_fmin())
    print(result)
    np.save('best_ARA_values.npy',m.values)
    return(m.values)


interp = LoadPropTimes() #load spicecore propagation times
cal_times = LoadCalTimes() #load calpulser propagation times
cal_delays = {(8, 25): 310.7580988219895, (16, 17): 73.63767997382199, (1, 24): 46.942653795811516, (24, 16): 0.3178174083769634, (25, 16): -128.64406086387436, (0, 25): 315.0504744764398, (1, 17): 120.86690117801047, (8, 9): 136.76464332460733, (9, 16): 45.349394633507856, (0, 16): 186.42203861256544, (1, 9): 2.004826570680628, (8, 17): 255.68921793193718, (0, 8): 4.308000654450262, (9, 25): 173.9465804973822, (25, 17): -55.084505890052355, (0, 1): 139.06781740837695, (24, 25): 128.9462532722513, (8, 24): 181.76497054973822, (1, 25): 175.93578206806282, (1, 16): 47.369846204188484, (0, 24): 186.0885962041885, (0, 17): 259.9815935863874, (24, 17): 73.87737238219896, (1, 8): -134.8066917539267, (9, 17): 118.90894960732984, (0, 9): 141.0882689790576, (9, 24): 44.98470222513089, (8, 16): 182.09841295811518}
cal_delays = {(0, 1): 138.54556609947645, (0, 8): 4.236174738219895, (0, 9): 140.63718913612564, (0, 24): 186.39258835078533, (0, 25): 315.067817408377, (0, 16): 186.48813808900525, (0, 17): 259.74550065445027, (1, 8): -134.35626636125656, (1, 9): 2.0759980366492146, (1, 24): 47.8001472513089, (1, 25): 176.47537630890054, (1, 16): 47.958196989528794, (1, 17): 121.18430955497382, (8, 9): 136.38538939790575, (8, 24): 182.14078861256544, (8, 25): 310.81601767015707, (8, 16): 182.26758835078533, (8, 17): 255.52495091623035, (9, 24): 45.73977421465968, (9, 25): 174.4150032722513, (9, 16): 45.866573952879584, (9, 17): 119.12393651832461, (24, 25): 128.65960405759162, (24, 16): 0.07992473821989529, (24, 17): 73.36853730366492, (25, 16): -128.59530431937173, (25, 17): -55.3379417539267, (16, 17): 73.30423756544502}
cal_delays = np.load('ARA_CalPulser_delaydict.npy',allow_pickle='TRUE').item()
pairs=list(itertools.combinations((0,1,8,9,24,25,16,17),2))
int_delays, SPICE_depths = LoadSPICEwf(pairs)



if __name__ == '__main__':
    first_vals = FirstIteration()
    second_vals = SecondIteration(first_vals)