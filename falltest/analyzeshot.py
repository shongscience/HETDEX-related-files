# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 10:25:37 2016

@author: shong
"""

import matplotlib.pyplot as plt
import numpy as np

import scipy.integrate as intg
from scipy.integrate import simps
from scipy.integrate import romberg
from scipy.interpolate import interp1d
from mpmath import gamma as Gamma
import scipy.optimize as op


""" ============================
Read ra,dec,shotname
"""
ras, decs, names = np.loadtxt("shot.dat",usecols=(1,2,8), unpack=True, dtype=np.str)

ra = ras.astype(np.float)
dec = decs.astype(np.float)

numshot = np.int(len(names))
name = np.empty(numshot,dtype='object')
track = np.empty(numshot,dtype='object')
trackflag = np.zeros(numshot)

for i in range(numshot):
    name[i] = names[i].split('_')[0]
    track[i] = names[i].split('_')[1]

iwest = np.where(track == 'W')
ieast = np.where(track == 'E')
trackflag[iwest] = 1


""" ============================
observed shot name
"""
obsdname = np.loadtxt("observedlist.dat",unpack=True, dtype=np.str)

obsdflag = np.zeros(numshot)

numobsdshot = len(obsdname)
iobsd = np.array([])

for i in range(numobsdshot):
    itmp = np.argwhere(name == obsdname[i])
    iobsd = np.append(iobsd,itmp)

iobsd = iobsd.astype(np.int)



""" ============================
plot : need ipython
"""
#plot basic settings
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman') 
plt.rcParams.update({'font.size': 16})

fig = plt.figure(figsize=(9,7))

plt.title('Fall Shots',fontsize=18)
plt.axis([1.335,1.405,-0.35,0.35])
plt.scatter(ra,dec,color='black',marker='.')
plt.scatter(ra[iobsd],dec[iobsd],s=50, color='red',marker='x')

offset = 0.001
for i in range(numshot):
    if i in iobsd:
        plt.text(ra[i]+offset,dec[i],name[i],color='red')
    else:
        plt.text(ra[i]+offset,dec[i],name[i])



#plt.axis('equal')
#plt.axes().set_aspect(1.0)
plt.xlabel(r'$R.A.$')
plt.ylabel(r'$Decl.$')


plt.savefig('fallshot.png')

plt.show()



