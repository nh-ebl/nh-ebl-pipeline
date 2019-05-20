#!/usr/bin/python

import numpy  as np
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    #return a*np.exp(-b*np.sqrt(x))+c#+d*np.sqrt(x)+e*x
    return a/x**2 + b
    
t,d = np.loadtxt('d_vs_t.txt',delimiter=',',unpack=True)

ds,tm = np.loadtxt('tm_vs_d.txt',delimiter=',',unpack=True)

mastert = np.arange(2006,2030,0.05)

t_transform = np.polyfit(t,d,1)

popt, pcov = curve_fit(func,ds,tm,p0=(1e6,2e1))

fig=plt.figure(figsize=(6.5,5))
ax = fig.add_subplot(1,1,1)

#ax.plot(t,d,linestyle=' ',marker='o',label='Temperature')
#ax.plot(mastert,t_transform[0] * mastert + t_transform[1])

#ax.plot(ds,tm,linestyle=' ',marker='o')
ax.semilogy(mastert,func(t_transform[0]*mastert+t_transform[1],popt[0],popt[1]),linestyle='-',label='Three-axis pointing Mode')
ax.semilogy(mastert,2.*func(t_transform[0]*mastert+t_transform[1],popt[0],popt[1]),linestyle='-',label='Spin-stabilized Mode')
ax.set_ylim([100,1e6])
ax.set_xlim([2006.1,2030])
ax.set_ylabel(r'Dual-TWTA DSN 70 m Data Rate (bps)',size=16)
ax.set_xlabel(r'Year',size=16)

ax.semilogy([2015.534,2015.534],[100,1e6],linestyle=':',color='k')
ax.annotate('Pluto Fly-by', xy=(2014.7,5e4), xytext=(2014.7,5e4),rotation=90)
            
ax.semilogy([2019.0,2019.0],[100,1e6],linestyle=':',color='k')
ax.annotate('2014 MU$_{69}$ Fly-by', xy=(2018.1,8e4), xytext=(2018.1,8e4),rotation=90)

ax.semilogy([2022,2022],[100,1e6],linestyle='--',color='red')
bbox_props = dict(boxstyle="rarrow,pad=0.3", fc="none", ec="red", lw=1)
t = ax.text(2025.57, 1e4, "Proposed Survey", ha="center", va="center",
            size=12,
            bbox=bbox_props)
#ax.annotate('Hypothetical Survey', xy=(2018.1,8e4), xytext=(2018.1,8e4))

plt.legend(loc=3,fontsize=12)

ax2 = ax.twiny()
ax2.semilogy(ds,tm,color='none')
ax2.set_ylim([100,1e6])
ax2.set_xlim([1.9,80])
ax2.set_xlabel(r'Heliocentric Distance (AU)',size=16)



#plt.tight_layout()
#plt.show()
plt.savefig('fig_datarate.pdf',bbox_inches='tight')

