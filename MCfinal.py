#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 15:58:13 2021

@author: aparnap

"""

import scipy as sp
import scipy.stats as stats
from random import random
import math

class MCsim():
    
    """
    This class contains three different Monte Carlo simulations for various different investigations
    the behaviour of gamma rays. Each needs a time input to determine the number of data points it will
    create given a certain activity. 
    """
    
    def __init__(self, time):
        self.t = time
        
    def single(self,energy,angle,activity):
        
        """
        This funciton distributes data points around a single energy value dependent on the 
        scattering angle using Monte Carlo methods and the Compton scattering formula. 
        The number of points within a certain energy range or "bin" are counted for use in data analysis.
        Input parameters are the unscattered gamma particle energy, the angle that the particles are
        scattered by and the activity of the source. Various constant parameters are contained within
        the function.
        """
        
        self.count = math.floor(activity*self.t)
        a=3355.55979142 # a coefficient for resolution vs energy relationship
        b=-16.72145688 # another coefficient for resolution vs energy relationship
        m = 511 # e- rest mass (keV)
        g = 5e-3 #detect gain (V per keV)
        bd = 9 #MCA bit depth
        maxs = 5.091 #MCA max signal (V)
        self.cspE = energy/(1+(energy/(m)*(1-sp.cos(angle*sp.pi/180)))) #scattered photon energy (keV)
        #print(cspE)
        res = sp.log(self.cspE/a)/b #calculation for resolution for scattered energy peak
        asig = g*self.cspE #analog signal (V) 
        #print(asig)
        rand_s=[] #empty random seed list
        for i in range(self.count):
            x = random()
            rand_s.append(x)#generates our set of random numbers and adds to list
        n = stats.norm.ppf(rand_s,0,res*asig/2) #distributes about ideal signal
        self.outb=[] #empty list to append final output energies into
        for i in range(len(rand_s)):
            x=math.floor(2**bd*(asig+n[i])/(maxs*0.5))
            self.outb.append(x)
        self.b=sp.arange(0,2001,5) #creating our energy bins
        self.count = sp.zeros(len(self.b))
        for i in range(len(self.b)):
            for j in range(len(self.outb)):
                if self.b[i] > self.outb[j] and self.b[i-1] <= self.outb[j]: #ascertains which bin to count a point in
                    self.count[i-1]+=1
        
    def scatter(self, energy, angle, activity):
        """
        Generates points around a scattered gamma energy peak and on a distribution representing
        the background data for a given gamma particle energy, scattering angle and source activity
        using Monte Carlo methods and Compton scattering. 90.7% of the data points are allocated to 
        the background and the remainder to the scattered peak. Again this funciton counts the number 
        of data points which occur in a certain energy range.
        """
        Emission_Peak = activity*self.t
        Background = 9.75*Emission_Peak
        self.num = math.floor(Emission_Peak + Background)     
        a=3355.55979142 # a coefficient for resolution vs energy relationship
        b=-16.72145688 # another coefficient for resolution vs energy relationship
        m = 511 # e- rest mass (keV)
        g = 5e-3 #detect gain (V per keV)
        bd = 9 #MCA bit depth
        maxs = 5.091 #MCA max signal (V)
        self.cspE = energy/(1+(energy/(m)*(1-sp.cos(angle*sp.pi/180)))) #scattered photon energy (keV)
        res = sp.log(self.cspE/a)/b #calculation for resolution at energy after scattering
        asig = g*self.cspE #analog signal (V) 
        rand_s=[] #empty random seed list
        for i in range(self.num):
            x = random()
            rand_s.append(x)#generates our set of random numbers and adds to list
        c=round(0.907*self.num)
        rb = rand_s[0:c] #random probabilities for background
        rc = rand_s[c:self.num] #randomly generated probabilities for scattered peak
        n = stats.norm.ppf(rc,0,res*asig/2) #distributes about 0 creating noise in peak
        outb=[] #empty list to append final output energies into
        for i in range(len(rc)):
            x=math.floor(2**bd*(asig+n[i])/maxs/0.5)
            outb.append(x)
        for i in range(len(rb)):
            outb.append(stats.gamma.ppf(rb[i],2,loc=0,scale=105)) #generating energies from background pdf
        self.b=sp.arange(0,1825,5) #creating our energy bins
        self.count = sp.zeros(len(self.b)) 
        for i in range(len(self.b)):
            for j in range(len(outb)):
                if self.b[i] > outb[j] and self.b[i-1] <= outb[j]: #ascertains which bin to count a point in
                    self.count[i-1]+=1
                    
    def peak2(self,energy1,energy2, activity):
        """
        Uses normal distributions to produce data around the two input energy values for a given activity. 
        This then counts the number of data points which occur within a given energy range. 
        """
        self.count = math.floor(2*activity*self.t) #assuming same activity and counts split 50/50
        a=3355.55979142 # a coefficient for resolution vs energy relationship
        b=-16.72145688 # another coefficient for resolution vs energy relationship
        g = 5e-3 #detect gain (V per keV)
        bd = 9 #MCA bit depth
        maxs = 5.091 #MCA max signal (V)
        res1= sp.log(energy1/a)/b #resolution at energy 1
        res2 = sp.log(energy2/a)/b #resolution at energy 2
        asig1 = g*energy1 #analog signal 1 (V) 
        asig2 = g*energy2 #analog signal 2 
        #print(res1, res2, E1, E2, asig1, asig2)
        rand_s=[] #empty random seed list
        for i in range(self.count):
            x = random()
            rand_s.append(x) #generates our set of random numbers and adds to list
        c=round(self.count/2)
        r1 = rand_s[0:c] #half list for first energy
        r2 = rand_s[c:self.count] # half list for second energy
        n1 = stats.norm.ppf(r1,0,res1*asig1/2) #distributes about ideal signal for energy 1
        n2 = stats.norm.ppf(r2,0,res2*asig2/2) #distributes about ideal signal for energy 1
        #print(n2)
        self.outb=[]
        for i in range(len(r1)):
            x=math.floor(2**bd*(asig1+n1[i])/(maxs*0.5))
            self.outb.append(x)
        for i in range(len(r2)):
            x=math.floor(2**bd*(asig2+n2[i])/(maxs*0.5))
            self.outb.append(x)
        #print(len(self.outb))
        self.b=sp.arange(0,2001,5) #generating bins of range 5keV
        self.count = sp.zeros(len(self.b))
        for i in range(len(self.b)):
            for j in range(len(self.outb)):
                if self.b[i] > self.outb[j] and self.b[i-1] <= self.outb[j]: #determining which bin to put data point in
                    self.count[i-1]+=1