#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 15:58:13 2021

@author: aparnap
"""

import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import MCfinal as mc

#%% Code used to collect data for scattering without background
E=[]
u=[]
s=mc.MCsim(1000) #input time
def g(x,a,b,c):
    y=a*sp.exp(-(x-b)**2/(2*c**2))
    return y
for i in range(5):
    s.single(662,20,1.7)
    b_=s.b[103:143] #for 20 degree scattering angle  
    c_=s.count[103:143] #for 20 degree scattering angle  
    #b_=s.b[93:133] #for 30 degree scattering angle 
    #c_=s.count[93:133] #for 30 degree scattering angle
    #b_=s.b[76:116] #for 45 degree scattering angle 
    #c_=s.count[76:116] #for 45 degree scattering angle 
    popt, pcov =curve_fit(g,b_,c_, [25,500,20])
    E.append(popt[1])
    u.append(pcov[1,1])
avE=sum(E)/5
avu=sum(u)/5
print(avE, avu)

#%% Preliminary data for scattered peak with background
s=mc.MCsim(2001) #input time
s.scatter(662,30,1.7)            
plt.scatter(s.b,s.count)
def g(x,a,b,c):
    y=a*sp.exp(-(x-b)**2/(2*c**2))
    return y
#These energy ranges are only valid for Cs-137
#b_=s.b[103:143] #for 20 degree scattering angle  
#c_=s.count[103:143] #for 20 degree scattering angle 
b_=s.b[93:133] #for 30 degree scattering angle 
c_=s.count[93:133] #for 30 degrees scattering angle 
#b_=s.b[76:116] #for 45 degree scattering angle 
#c_=s.count[76:116] #for 45 degree scattering angle
popt, pcov =curve_fit(g,b_,c_, [25,500,20])
plt.plot(b_, g(b_, *popt),c='k')
print('simulation = ', popt[1], '+/-', pcov[1,1], ',', 'ideal = ', s.cspE)
plt.ylabel('Number of counts')
plt.xlabel('Energy (keV)')
#plt.xticks(fontsize=12.5, rotation=0)
#plt.yticks(fontsize=12.5, rotation=0)
plt.grid()
plt.show()
#%% Code used to take data for scattering with background
E=[]
u=[]
def g(x,a,b,c):
    y=a*sp.exp(-(x-b)**2/(2*c**2))
    return y
s=mc.MCsim(4256)
for i in range(5):
    s.scatter(662,45, 1.7)
    #b_=s.b[103:143] #for 20 degree scattering angle  
    #c_=s.count[103:143] #for 20 degree scattering angle  
    #b_=s.b[93:133] #for 30 degree scattering angle 
    #c_=s.count[93:133] #for 30 degree scattering angle
    b_=s.b[76:116] #for 45 degree scattering angle 
    c_=s.count[76:116] #for 45 degree scattering angle 
    popt, pcov =curve_fit(g,b_,c_, [25,500,20])
    E.append(popt[1])
    u.append(pcov[1,1])
avE=sum(E)/5
avu=sum(u)/5
print(avE, avu)
#%% 2 peak detection testing
s=mc.MCsim(489)
s.peak2(450,550,1.7)
plt.scatter(s.b,s.count)
plt.xlim(0,1000)
def g(x,a,b,c):
    y=a*sp.exp(-(x-b)**2/(2*c**2))
    return y
opt1, cov1 =curve_fit(g,s.b[50:100],s.count[50:100], [25,450,20])
plt.plot(s.b[50:100], g(s.b[50:100], *opt1),c='k')
opt2, cov2 = curve_fit(g,s.b[100:150],s.count[100:150], [25,450,20])
plt.plot(s.b[100:150], g(s.b[100:150], *opt2),c='k')
plt.show()
print(opt1[1],cov1[1,1],opt2[1],cov2[1,1])
#%% Taking final data for 2 peaks
E1=[]
u1=[]
E2=[]
u2=[]
def g(x,a,b,c):
    y=a*sp.exp(-(x-b)**2/(2*c**2))
    return y
s=mc.MCsim(647)
for i in range(5):
    s.peak2(450,550,1.7)
    opt1, cov1 =curve_fit(g,s.b[50:100],s.count[50:100], [25,450,20])
    opt2, cov2 = curve_fit(g,s.b[100:150],s.count[100:150], [25,550,20])   
    E1.append(opt1[1])
    u1.append(cov1[1,1])
    E2.append(opt2[1])
    u2.append(cov2[1,1])
avE1=sum(E1)/5
avu1=sum(u1)/5
avE2=sum(E2)/5
avu2=sum(u2)/5
print(avE1, avu1,avE2, avu2)
