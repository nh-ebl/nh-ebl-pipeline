
# coding: utf-8

# In[2]:

get_ipython().magic(u'matplotlib inline')


# In[3]:

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy

from os import path
from sys import argv


# ### DEFINE FUNCTIONS

# In[4]:

def load_data(data_filename):
    # load data file
    fields = np.genfromtxt(data_filename, delimiter=',',unpack=True,usecols=(0))
    ID = np.genfromtxt(data_filename, delimiter=',',unpack=True,usecols=(1),dtype=np.str)
    lorri_m = np.genfromtxt(data_filename, delimiter=',',unpack=True,usecols=(2))
    S, S_err = np.genfromtxt(data_filename, delimiter=',',unpack=True,usecols=(3,4))
    m, m_err = np.genfromtxt(data_filename, delimiter=',',unpack=True,usecols=(5,6))
    
    return (fields,ID,lorri_m,S,S_err,m,m_err)

def find_weighted_average(data,err):
    weight = 1.0/(err**2.0) # weight the data points by their errors. Smaller error = more weight
    weight = weight/np.sum(weight) # normalize the weight
    
    mean = np.average(data, weights=weight) 
    #calculate the unbiased variance of the mean
    mean_var_1 = np.sum(weight*(data-mean)**2.0)
    mean_var_2 = np.sum(weight)/((np.sum(weight))**2.0-np.sum(weight*weight))
    mean_var = mean_var_1*mean_var_2
    #error = sqrt of the variance
    mean_err = mean_var**(1.0/2)
    
    return (mean,mean_err)

def find_average(data): 
    mean = np.average(data)
    err = np.std(data)
    
    return (mean,err)
    
def find_flux(lorri_m,S,S_err):
    a = 10**(-lorri_m/2.5)
    flux = a / S # unit: DN
    flux_err = np.sqrt(a*a*S_err*S_err/(S**4.0)) # unit: DN
    
    #convert DN to Jy in R_L band
    flux = flux*3050
    flux_err = flux_err*3050
    
    #return flux in microJy
    flux = flux*(10**6.0)
    flux_err = flux_err*(10**6.0)
    
    return (flux,flux_err)


# ### LOAD DATA

# In[5]:

all_fields_data = 'ALL_FIELDS.csv'
good_fields_data = 'GOOD_FIELDS.csv'

ALL_fields,ALL_ID,ALL_lorri_m,ALL_S,ALL_S_err,ALL_m,ALL_m_err = load_data(all_fields_data)
GOOD_fields,GOOD_ID,GOOD_lorri_m,GOOD_S,GOOD_S_err,GOOD_m,GOOD_m_err = load_data(good_fields_data)

#Values from Cheng et al. to compare our numbers to
CHENG_m = 18.86 
CHENG_flux = 10**(-CHENG_m/2.5)*3050*(10**(6.0)) #assume a count of 1 per sec, then mag measured = ref mag = 18.86
                                     #3050 is to convert to Jy, 10^6 is to convert to microJy
print('Cheng et al. mag: %0.3f' % CHENG_m)
print('Cheng et al. flux: %0.3f [microJy]' % CHENG_flux)


# ### COMPUTE MEAN MAG AND FLUX

# In[6]:

#MAGNITUDE
#only good fields
GOOD_m_mean, GOOD_m_mean_err = find_weighted_average(GOOD_m,GOOD_m_err)
print("Mean mag of 4 best fields: %0.3f +/- %0.3f" % (GOOD_m_mean,GOOD_m_mean_err))

#all fields
ALL_m_mean, ALL_m_mean_err = find_average(ALL_m)
print("Mean mag of all fields: %0.3f +/- %0.3f" % (ALL_m_mean,ALL_m_mean_err))


# In[7]:

#FLUX
# a = 10^(-lorri_m/2.5)
# F = 3050*a/s
# <F> = <a>/<s>, <a> = <10^(-lorri_m/2.5)> ~ = 10^(-<lorri_m>/2.5)
# <lorri_m> is the average of lorri_m of EACH star, with associated stdev
# so <lorri_m> = lorri_m reported here
# <S> = S reported here, so the flux of each star will then be:
# <F_each star> = 3050*a/S
# error of each point is (use derivative to get the average):
# <F_err_each star> = sqrt(a_err*a_err + a*a/(S_err*S_err)*S_err)/S

# THEN use this F_err to find the weighted mean

#only good fields
GOOD_flux,GOOD_flux_err = find_flux(GOOD_lorri_m,GOOD_S,GOOD_S_err)
GOOD_flux_mean,GOOD_flux_mean_err = find_weighted_average(GOOD_flux,GOOD_flux_err)
print("Mean flux of 4 best fields: %0.3f +/- %0.3f [microJy]" % (GOOD_flux_mean,GOOD_flux_mean_err))

#all fields
ALL_flux,ALL_flux_err = find_flux(ALL_lorri_m,ALL_S,ALL_S_err)
ALL_flux_mean,ALL_flux_mean_err = find_average(ALL_flux)
print("Mean flux of all fields: %0.3f +/- %0.3f [microJy]" % (ALL_flux_mean,ALL_flux_mean_err))


# ### MAKING FIGURES

# In[8]:

#EXTRACT THE POINTS OF FOUR BEST FIELDS FOR PLOT
FIELD_3 = np.asarray([1,2])
FIELD_5 = np.asarray([3,4])
FIELD_6 = np.asarray([5,6])
FIELD_7 = np.asarray([7,8])
star = np.arange(1,9,1)

FIELD_3_m = GOOD_m[0:2]
FIELD_5_m = GOOD_m[2:4]
FIELD_6_m = GOOD_m[4:6]
FIELD_7_m = GOOD_m[6:]

FIELD_3_m_err = GOOD_m_err[0:2]
FIELD_5_m_err = GOOD_m_err[2:4]
FIELD_6_m_err = GOOD_m_err[4:6]
FIELD_7_m_err = GOOD_m_err[6:]

FIELD_3_flux = GOOD_flux[0:2]
FIELD_5_flux = GOOD_flux[2:4]
FIELD_6_flux = GOOD_flux[4:6]
FIELD_7_flux = GOOD_flux[6:]

FIELD_3_flux_err = GOOD_flux_err[0:2]
FIELD_5_flux_err = GOOD_flux_err[2:4]
FIELD_6_flux_err = GOOD_flux_err[4:6]
FIELD_7_flux_err = GOOD_flux_err[6:]


# In[10]:

#PLOT
fig,(ax1,ax2) = plt.subplots(1,2,figsize=(9,4.5))

#MAGNITUDE
star_number = np.arange(0,10,1) # x-coordinate for some points to use fill_between

#Cheng et al.
ax1.axhline(CHENG_m,color='#fe9929',linestyle='-',linewidth=0.75,zorder=0,label='Cheng et al.')

#all stars average
ax1.axhline(ALL_m_mean,color='#1d91c0',linestyle=':',linewidth=1,zorder=1,label='All fields')

#four fields average
ax1.axhline(GOOD_m_mean,color='#41ae76',linestyle='-',linewidth=2,label='4 best fields',zorder=2,alpha=0.7)
ax1.fill_between(star_number,GOOD_m_mean+GOOD_m_mean_err,GOOD_m_mean-GOOD_m_mean_err,                color='#66c2a4',alpha=0.15,zorder=2)

#individual stars
ax1.errorbar(FIELD_3,FIELD_3_m,yerr=FIELD_3_m_err,fmt='o',color='#c51b7d',ms=4,zorder=6,linewidth=0.7)
ax1.errorbar(FIELD_5,FIELD_5_m,yerr=FIELD_5_m_err,fmt='^',color='#c51b7d',ms=5,zorder=6,linewidth=0.7)
ax1.errorbar(FIELD_6,FIELD_6_m,yerr=FIELD_6_m_err,fmt='s',color='#c51b7d',ms=4,zorder=6,linewidth=0.7)
ax1.errorbar(FIELD_7,FIELD_7_m,yerr=FIELD_7_m_err,fmt='*',color='#c51b7d',ms=4,zorder=6,linewidth=0.7)

ax1.set_xlabel('Star ID',labelpad=8)
ax1.set_ylabel(r'Reference magnitude $R_{{\rm L},0}$')
ax1.set_xlim(0,9)
ax1.set_ylim(18.1,19.1)
ax1.set_xticks(star) # set the ticks to star number
ax1.set_xticklabels(GOOD_ID,rotation=90) # match star to star ID
ax1.legend(loc='lower left',frameon=False, prop={'size':10})


#FLUX
F_mean_star_plot = np.arange(0,9,0.5)

#Cheng et al.
ax2.axhline(CHENG_flux,color='#fe9929',linestyle='-',linewidth=0.75,zorder=0,label='Cheng et al.')

#all stars average
ax2.axhline(ALL_flux_mean,color='#1d91c0',linestyle=':',linewidth=1,zorder=1,label='All fields')

#four fields average
ax2.axhline(GOOD_flux_mean,color='#41ae76',linestyle='-',linewidth=2,label='4 best fields',zorder=2,alpha=0.7)
ax2.fill_between(star_number,GOOD_flux_mean+GOOD_flux_mean_err,GOOD_flux_mean-GOOD_flux_mean_err,                 color='#66c2a4',alpha=0.15,zorder=2)

#individual stars
ax2.errorbar(FIELD_3,FIELD_3_flux,yerr=FIELD_3_flux_err,fmt='o',color='#c51b7d',ms=4,zorder=6,linewidth=0.7)
ax2.errorbar(FIELD_5,FIELD_5_flux,yerr=FIELD_5_flux_err,fmt='^',color='#c51b7d',ms=5,zorder=6,linewidth=0.7)
ax2.errorbar(FIELD_6,FIELD_6_flux,yerr=FIELD_6_flux_err,fmt='s',color='#c51b7d',ms=4,zorder=6,linewidth=0.7)
ax2.errorbar(FIELD_7,FIELD_7_flux,yerr=FIELD_7_flux_err,fmt='*',color='#c51b7d',ms=4,zorder=6,linewidth=0.7)

ax2.set_xlabel('Star ID',labelpad=8)
ax2.set_ylabel(r'Reference flux $F_{{\rm L}, 0}$ [$\mu$Jy]')
ax2.set_xlim(0,9)
ax2.set_ylim(55,155)
ax2.set_xticks(star)
ax2.set_xticklabels(GOOD_ID,rotation=90)
ax2.legend(loc='lower right',frameon=False, prop={'size':10})



plt.subplots_adjust(bottom=0.3)
plt.tight_layout(w_pad=5.0)
plt.savefig('aperture_calibration.pdf')
plt.show()

