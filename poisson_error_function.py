#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:41:36 2017

@author: eejvt
"""

import numpy as np
from random import random
from random import uniform
import random
import numpy.random as nprnd
#import Jesuslib as jl
import matplotlib.pyplot as plt
from glob import glob

def find_nearest_vector_index(array, value):
    n = np.array([abs(i-value) for i in array])
    nindex=np.apply_along_axis(np.argmin,0,n) 

    return nindex



def calc_ff(INP_expected_value,droplets):
    #drop_array=np.zeros(droplets)
    INPs=np.random.poisson(INP_expected_value,1)
    pos_inp=nprnd.randint(droplets,size=INPs)#position of every INP in droplets
    frezzed_droplets=0
    frezzed_droplets=len(set(pos_inp.tolist()))
    return float(frezzed_droplets)/float(droplets)


    
def poisson_errors(N, mu=1,sigma=0):
    Nislist=isinstance(N, list)
    
    
    #def errors(N):
        
    
#    file_name='/nfs/see-fs-01_users/eejvt/Run_labradorite.csv'
#    data=np.genfromtxt(file_name,delimiter=',',skip_header=1)
#    ts=data[:,0]
#    print file_name
    
    
    
    if Nislist:ts=N
    else:ts=np.arange(N)
    
    if mu==1 and sigma==0 and N<100:
        try:
            final_data=np.genfromtxt('/nfs/see-fs-01_users/eejvt/PYTHON_CODE/INP_poisson_errors/pre_calculated/%i.csv'%N,skip_header=1,delimiter=',')
            print 'Poisson errors read from precalculated values. \n Droplets:',N
            return final_data
        except:
            print 'simulating poisson'
        
    #mu, sigma = data[0,1], data[0,2]/1.96# mean and standard deviation#1.96
    
    droplets=len(ts)
    
    #droplets=24
    #mu=1
    #sigma=0
    
    
    print 'Droplets %i'%droplets
    plot=0
    if plot:
        plt.figure(figsize=(15,12))
    target_ff=np.linspace(0,1,droplets+1)[1:]
    repeat=1000
    steps=2500
    particles_per_droplet=np.linspace(1e-2,5,steps)
    particles_per_droplet=np.logspace(-4,1.3,steps)
    
    particles_in_total_array=particles_per_droplet*droplets
    ff_arr=np.zeros((repeat,particles_in_total_array.size))
    mean_ff=np.zeros(particles_in_total_array.size)
    lower_95=np.zeros(particles_in_total_array.size)
    upper_95=np.zeros(particles_in_total_array.size)
    mean_ff=np.zeros(target_ff.size)
    lower_95=np.zeros(target_ff.size)
    upper_95=np.zeros(target_ff.size)
    
    
    for iparticles in range(particles_in_total_array.size):
        #print iparticles
        for i in range(repeat):
            ff_arr[i,iparticles]=calc_ff(particles_in_total_array[iparticles],droplets)
    
    
    all_ff_list=np.zeros(target_ff.size).tolist()
    cumulative=np.zeros((steps,target_ff.size))
    print 'poisson simulated'
    for iff in range(target_ff.size):
        ff_list=[]
        for istep in range(steps):
            for iv in range(np.array([np.abs(ff_arr[:,istep]-target_ff[iff])<0.0001]).sum()):
                ff_list.append(particles_per_droplet[istep])
            cumulative[istep,iff]=np.array([np.abs(ff_arr[:,istep]-target_ff[iff])<0.0001]).sum()
        #mean_ff[iff]=np.array(ff_list).mean()
        mean_ff[iff]=np.mean(ff_list)
        mean_ff[iff]=np.sort(ff_list)[int(len(ff_list)*0.5)]
        if len(ff_list)!=0:
            lower_95[iff]=np.sort(ff_list)[int(len(ff_list)*0.025)]
            upper_95[iff]=np.sort(ff_list)[int(len(ff_list)*0.975)]
            target_ff[iff]
        all_ff_list[iff]=np.array(ff_list)
    final_data=np.vstack([target_ff,mean_ff,lower_95,upper_95]).T
    cumulative_normalized=cumulative/(repeat*steps)
    
    
    
    
    cumulative_ff_normalized=np.zeros(cumulative.shape)
    for i in range(cumulative.shape[1]):
        cumulative_ff_normalized[:,i]=cumulative[:,i]/float(np.sum(cumulative[:,i]))
    import matplotlib.pyplot as plt
    from scipy import stats
    
    print 'propagating errors'
    montecarlo_values=10000
    final_data=np.zeros((len(target_ff),5))
    final_data[:,0]=target_ff
    final_data[:,1]=ts
    
    for ff_index in range(len(target_ff)):
        print ff_index+1,'out of', len(target_ff)
        xk = np.arange(0,len(particles_per_droplet),1)
        pk = cumulative_ff_normalized[:,ff_index]
        custm = stats.rv_discrete(name='custm', values=(xk,pk))
        #h = plt.plot(xk, custm.pmf(xk))
        R = custm.rvs(size=montecarlo_values)
        k_values=particles_per_droplet[R]
        if sigma>0:
            factor_values = np.random.normal(mu, sigma, montecarlo_values)
        else:
            factor_values=np.ones(montecarlo_values)*mu
        ns_values=k_values/factor_values
        ns_lower_95=np.sort(ns_values)[int(len(ns_values)*0.025)]
        ns_upper_95=np.sort(ns_values)[int(len(ns_values)*0.975)]
        ns_median=np.sort(ns_values)[int(len(ns_values)*0.5)]
        
        ns_mean=-np.log(1-target_ff[ff_index])/mu
        #plt.hist(ns_values, 3000, normed=True)
        #plt.axvline(ns_lower_95,color='b',lw=2)
        #plt.axvline(ns_upper_95,color='b',lw=2)
        #plt.axvline(ns_mean,color='r',lw=2)
        #plt.axvline(ns_median,color='g',lw=2)
        final_data[ff_index,2]=ns_mean
        final_data[ff_index,3]=ns_upper_95
        final_data[ff_index,4]=ns_lower_95
    

    return final_data
    


