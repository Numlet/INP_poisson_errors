#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 10:01:58 2017

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


#def calc_ff(INP,droplets):
#    #drop_array=np.zeros(droplets)
#    if (INP-int(INP))>random.random():
#        INP=float(int(INP)+1)
#    else:
#        INP=float(int(INP))
#    pos_inp=nprnd.randint(droplets,size=INP)#position of every INP in droplets
#    frezzed_droplets=0
#    frezzed_droplets=len(set(pos_inp.tolist()))
#    return float(frezzed_droplets)/float(droplets)


N=24
mu=1
sigma=0


Nislist=(len([N])>1)


#def errors(N):
    


file_name='/nfs/see-fs-01_users/eejvt/Run_labradorite.csv'
data=np.genfromtxt(file_name,delimiter=',',skip_header=1)
ts=data[:,0]
print file_name



if Nislist:ts=N
else:ts=np.arange(N)

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
    
    
np.savetxt('/nfs/see-fs-01_users/eejvt/poisson_test_%i_droplets.csv'%droplets,final_data,fmt='%.6f',delimiter=',',header='Fraction frozen, Temperature, ns,ns_upper,ns_lower, A=%.6f,std=%.6f'%(mu,sigma))
print 'final file saved'


plt.figure()
teoric_k=-np.log(1-target_ff)
teoric_k[-1]=teoric_k[-2]
plt.plot(lower_95,target_ff,'yo',lw=2)
plt.plot(upper_95,target_ff,'bo',lw=2)
plt.xscale('log')
#plt.yscale('log')
plt.errorbar(mean_ff,target_ff,xerr=[mean_ff-lower_95,upper_95-mean_ff])
plt.plot(mean_ff,target_ff,'g-',lw=2)
plt.ylim(0,1.2)#
plt.title('%i Droplets'%droplets)
plt.xlabel('K exponent')
plt.ylabel('Fraction Frozen')
plt.plot(teoric_k,target_ff,'r-',lw=2,label='teorical ff')
#return 



#%%




plt.figure()
bins=np.linspace(1000,250000,250)
plt.hist(np.random.poisson(100000,10000),bins=bins,alpha=.3,normed=True)

propagatedN=[]
propagated9=[]
propagated7=[]
propagated5=[]
propagated3=[]
propagated2=[]
propagated1=[]
for i in range(10000):
    print i
    mean1=np.random.poisson(100000,1)
    mean2=np.random.poisson(mean1,1)
    mean3=np.random.poisson(mean2,1)
    mean4=np.random.poisson(mean3,1)
    mean5=np.random.poisson(mean4,1)
    mean6=np.random.poisson(mean5,1)
    mean7=np.random.poisson(mean6,1)
    mean8=np.random.poisson(mean7,1)
    mean9=np.random.poisson(mean8,1)
    meanN=np.random.poisson(mean9,1)
    for _ in range(100):
        meanN=np.random.poisson(meanN,1)
    propagatedN.append(meanN)
    propagated9.append(mean9)
    propagated7.append(mean7)
    propagated5.append(mean5)
    propagated3.append(mean3)
    propagated2.append(mean2)
    propagated1.append(mean1)

plt.hist(np.array(propagatedN),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated9),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated7),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated5),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated3),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated2),bins=bins,alpha=.3,normed=True)
#plt.hist(np.array(propagated1),bins=np.linspace(1,250,250),alpha=.3,normed=True)

#%%
plt.figure()
bins=np.linspace(80,120,40)
plt.hist(np.random.normal(100,size=10000),bins=bins,alpha=.3,normed=True)

propagated9=[]
propagated7=[]
propagated5=[]
propagated3=[]
propagated2=[]
propagated1=[]
for i in range(10000):
    print i
    mean1=np.random.normal(100)
    mean2=np.random.normal(mean1)
    mean3=np.random.normal(mean2)
    mean4=np.random.normal(mean3)
    mean5=np.random.normal(mean4)
    mean6=np.random.normal(mean5)
    mean7=np.random.normal(mean6)
    mean8=np.random.normal(mean7)
    mean9=np.random.normal(mean8)
    propagated9.append(mean9)
    propagated7.append(mean7)
    propagated5.append(mean5)
    propagated3.append(mean3)
    propagated2.append(mean2)
    propagated1.append(mean1)

plt.hist(np.array(propagated9),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated7),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated5),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated3),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated2),bins=bins,alpha=.3,normed=True)
plt.hist(np.array(propagated1),bins=bins,alpha=.3,normed=True)
#plt.plot(ts,final_data[:,3],'k--',lw=2)
#plt.plot(ts,final_data[:,4],'k--',lw=2)
#plt.plot(ts,final_data[:,2],'k-',lw=2)
#plt.yscale('log')
##plt.yscale('log')
#plt.ylim(0,1.2)#
#plt.title('%i Droplets'%droplets)
#plt.xlabel('K exponent')
#plt.ylabel('Fraction Frozen')
#plt.plot(teoric_k,target_ff,'r-',lw=2)
#%%
particles=1000
boxes=1
arange=np.zeros(boxes)






#%%
#plt.errorbar(final_data[:,1],final_data[:,2],yerr=[final_data[:,2]-final_data[:,4],final_data[:,2]-final_data[:,3]])
'''
ns_values=final_data[:,2]
ns_up=final_data[:,3]
ns_low=final_data[:,4]
plt.errorbar(ts,ns_values,yerr=[ns_values-ns_low,ns_up-ns_values],c='k')
plt.yscale('log')
plt.title('%i Droplets'%droplets)
plt.xlabel('T')
plt.ylabel('$n_s$')
#%%
plt.plot(ts,final_data[:,3],'k--',lw=2)
plt.plot(ts,final_data[:,4],'k--',lw=2)
plt.plot(ts,final_data[:,2],'k-',lw=2)
plt.yscale('log')
#plt.yscale('log')
plt.ylim(0,1.2)#
plt.title('%i Droplets'%droplets)
plt.xlabel('K exponent')
plt.ylabel('Fraction Frozen')
plt.plot(teoric_k,target_ff,'r-',lw=2)
'''
#%%
'''
import matplotlib.pyplot as plt
count, bins, ignored = plt.hist(s, 3000, normed=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='r')
'''
#%%

'''
fig=plt.figure()
#ff_index=6
ax=plt.subplot(111)
plt.plot(particles_per_droplet,cumulative_normalized[:,ff_index],'ko',lw=0.5)
plt.plot(particles_per_droplet,cumulative_normalized[:,ff_index],'k--')
plt.axvline(lower_95[ff_index],color='r',lw=2)
plt.axvline(upper_95[ff_index],color='r',lw=2)
plt.axvline(mean_ff[ff_index],color='g',lw=2)
plt.xscale('log')
for r in R:
    plt.axvline(particles_per_droplet[r],color='b',lw=0.03)
'''




