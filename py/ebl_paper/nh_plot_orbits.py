"""
Demo of a line plot on a polar axis.
"""
import numpy as np
import matplotlib.pyplot as plt
#import png
#import itertools
from scipy.io import loadmat
#from scipy.ndimage import filters,zoom
#import skimage
#from skimage import io,color
#from PIL import Image
#import pyfits
from matplotlib.backends.backend_pdf import PdfPages

oneau = 1.495978706137e8

# set up a plot

mytraj = np.recfromcsv('lookup/nh_mission_trajectory.tab.txt',delimiter=',')

myearth = np.recfromcsv('lookup/earth_2006.txt',delimiter=',')
mymars = np.recfromcsv('lookup/mars_2006.txt',delimiter=',')
myjupiter = np.recfromcsv('lookup/jupiter_2007.txt',delimiter=',')
mysaturn = np.recfromcsv('lookup/saturn_1980.txt',delimiter=',')
myanus = np.recfromcsv('lookup/uranus_1950.txt',delimiter=',')

nhpos = loadmat('../../scratch/nh_position.mat')

mytraj['y'] = mytraj['y'] + 0.075 * oneau

mytrajr = np.sqrt(mytraj['x']**2 + mytraj['y']**2)
myearthr = np.sqrt(myearth['x']**2 + myearth['y']**2)
mymarsr = np.sqrt(mymars['x']**2 + mymars['y']**2)
myjupiterr = np.sqrt(myjupiter['x']**2 + myjupiter['y']**2)
mysaturnr = np.sqrt(mysaturn['x']**2 + mysaturn['y']**2)
myanusr = np.sqrt(myanus['x']**2 + myanus['y']**2)


fig=plt.figure(figsize=(10.9,5))

####################################################################

ax = fig.add_subplot(1,2,1)

ax.plot([0],[0],marker='o',color='yellow')

plt.tick_params(axis='both', which='major', labelsize=16)

ax.set_ylim([-22,22])
ax.set_xlim([-22,22])
ax.set_xlabel('$d_{x}$ (AU, Solar Ecliptic J2000.0)',fontsize=18)
ax.set_ylabel('$d_{y}$ (AU, Solar Ecliptic J2000.0)',fontsize=18)

#xm = np.squeeze(-nhpos.get('x','err'))
#ym = np.squeeze(nhpos.get('t','err'))
#whplnan = ~np.isnan(np.squeeze(nhpos.get('x','err')))
trajtime = (mytraj['time']-mytraj['time'][0])/(24*3600)

whdark = np.where(trajtime < 221.2)
ax.plot(mytraj['x'][np.squeeze(whdark)]/oneau,mytraj['y'][np.squeeze(whdark)]/oneau,linestyle='-',color='black',linewidth=3)

ax.plot(mytraj['x']/oneau,mytraj['y']/oneau,linestyle='-',color='red',linewidth=1.5)

whpl = np.abs(myjupiter['jdtdb'] - 2454155.500000000).argmin()

lfs = 13

ax.plot(myearth['x'],myearth['y'],linestyle=':',color='blue')
ax.plot(myearth['x'][19],myearth['y'][19],marker='x',color='black',markeredgewidth=2)
ax.plot(mymars['x'],mymars['y'],linestyle=':',color='blue')
ax.plot(myjupiter['x'],myjupiter['y'],linestyle=':',color='blue')
ax.plot(myjupiter['x'][whpl],myjupiter['y'][whpl],marker='x',color='black',markeredgewidth=2)
ax.plot(mysaturn['x'][0:140],mysaturn['y'][0:140],linestyle=':',color='blue')
ax.plot(myanus['x'][0:86],myanus['y'][0:86],linestyle=':',color='blue')

ax.annotate('Launch',xy=(-1,1),xycoords='data',
            xytext=(-5,5), textcoords='data',fontsize=lfs,
            arrowprops=dict(facecolor='black',arrowstyle='->'),
            horizontalalignment='right', verticalalignment='top',
            )

ax.annotate('Jupiter Encounter',xy=(-2.5,-4.8),xycoords='data',
            xytext=(-3,-1), textcoords='data',fontsize=lfs,
            arrowprops=dict(facecolor='black',arrowstyle='->'),
            horizontalalignment='right', verticalalignment='top',
            )

ax.annotate('Cover Ejected',xy=(-2,-2.5),xycoords='data',
            xytext=(3.5,0.5), textcoords='data',fontsize=lfs,
            arrowprops=dict(facecolor='black',arrowstyle='->'),
            horizontalalignment='left', verticalalignment='top',
            )

ax.annotate('Uranus',xy=(14,15),xycoords='data',\
             xytext=(14,15),textcoords='data',\
             color='blue',size=lfs)
ax.annotate('Saturn',xy=(7,7),xycoords='data',\
             xytext=(7,7),textcoords='data',\
             color='blue',size=lfs)
ax.annotate('Jupiter',xy=(4,4),xycoords='data',\
             xytext=(4,4),textcoords='data',\
             color='blue',size=lfs)
             
leaderlen = 8
ifield = 1

for i in range(0,np.size(nhpos['fieldflag'])):
    if nhpos['fieldflag'][i,0] == 1:# & ifield < 2 :
        print str(nhpos['fieldflag'][i,0])
        print str(nhpos['se'][i]) + ':' + str(nhpos['elong'][i])
        ax.arrow(np.squeeze(-nhpos['x'][i]),\
                  np.squeeze(-nhpos['y'][i]),\
                  np.squeeze(leaderlen * np.cos(nhpos['elong'][i]*np.pi/180)),\
                  np.squeeze(leaderlen * np.sin(nhpos['elong'][i]*np.pi/180)),\
                  linestyle='-',color='red',head_width=1, head_length=1,\
                  label=str(ifield))
        if np.squeeze(-nhpos['x'][i])+\
           1.0*np.squeeze(leaderlen * \
                          np.cos(nhpos['elong'][i]*np.pi/180)) < 0:
            padderx = 2.1
            paddery = 2.0
        else:
            padderx = 1.25
            paddery = 1.4
        
        ax.annotate('Field '+ str(ifield),\
                     xy=(np.squeeze(-nhpos['x'][i]),\
                         np.squeeze(-nhpos['y'][i])),\
                     xycoords='data',\
                     xytext=(np.squeeze(-nhpos['x'][i])+\
                             padderx*np.squeeze(leaderlen * \
                                        np.cos(nhpos['elong'][i]*np.pi/180)),\
                             np.squeeze(-nhpos['y'][i])+\
                             paddery*np.squeeze(leaderlen * \
                                        np.sin(nhpos['elong'][i]*np.pi/180))),\
                     textcoords='data',\
                     horizontalalignment='left',verticalalignment='center',\
                    fontsize=lfs)
                             
        ax.plot(-nhpos['x'][i],-nhpos['y'][i],color='red',marker='o')
        ifield = ifield + 1
        
####################################################################

ax1 = fig.add_subplot(1,2,2)

ax1.plot([0],[0],marker='o',color='yellow')

plt.tick_params(axis='both', which='major', labelsize=16)

ax1.set_ylim([-1,1])
ax1.set_xlim([0,26])

ax1.set_xlabel('$d_{r}$ (AU, Solar Ecliptic J2000.0)',fontsize=18)
ax1.set_ylabel('$d_{z}$ (AU, Solar Ecliptic J2000.0)',fontsize=18)

ax1.plot([0,26],[0,0],linestyle=':',color='black')



#nhpos['z'] = (nhpos['z']-0.32)/36

ax1.plot([myearthr[19],myearthr[19]],[-1,0],linestyle=':',color='blue')
ax1.plot([np.mean(mymarsr),np.mean(mymarsr)],[-1,0],linestyle=':',color='blue')
#ax.plot(myjupiter['x'],myjupiter['y'],linestyle=':',color='blue')
ax1.plot([myjupiterr[whpl],myjupiterr[whpl]],[-1,0.04],linestyle=':',color='blue')
ax1.plot([np.mean(mysaturnr[0:140]),np.mean(mysaturnr[0:140])],[-1,0.12],linestyle=':',color='blue')
ax1.plot([np.mean(myanusr[0:86]),np.mean(myanusr[0:86])],[-1,0.22],linestyle=':',color='blue')

ax1.annotate('Uranus',xy=(19.5,-0.7),xycoords='data',\
             xytext=(19.5,-0.7),textcoords='data',\
             color='blue',size=lfs,rotation='vertical')
ax1.annotate('Saturn',xy=(10,-0.7),xycoords='data',\
             xytext=(10,-0.7),textcoords='data',\
             color='blue',size=lfs,rotation='vertical')
ax1.annotate('Jupiter',xy=(5.6,-0.7),xycoords='data',\
             xytext=(5.6,-0.7),textcoords='data',\
             color='blue',size=lfs,rotation='vertical')

ax1.plot(mytrajr[np.squeeze(whdark)]/oneau,-(mytraj['z'][np.squeeze(whdark)]/oneau-0.32)/36,linestyle='-',color='black',linewidth=3)
ax1.plot(mytrajr/oneau,-(mytraj['z']/oneau-0.32)/36,linestyle='-',color='red',linewidth=1.5)
#plt.plot(mytrajr/oneau,mytraj['z']/oneau,marker='+')
#plt.plot(myearthr,myearth['z'],linestyle='-')
ax1.plot(myearthr[19],myearth['z'][19],marker='x',color='black',markeredgewidth=2)
#plt.plot([0,myjupiterr[whpl]],[0,myjupiter['z'][whpl]],linestyle=':')
ax1.plot(myjupiterr[whpl],myjupiter['z'][whpl],marker='x',color='black',markeredgewidth=2)

ax2 = ax1.twinx()
ax2.set_ylim([-13,13])
ax2.set_xlim([0,26])
ax2.get_yaxis().set_ticks([])
#ax2.set_ylabel('$b_{e}$ (degrees)')
#ax2.plot(mytrajr/oneau,np.sin(mytrajr/oneau))

#ax2.set_yticklabels([' ','-30','-15','0','15','30'])

ax1.annotate('Launch',xy=(0.98,0.02),xycoords='data',
             xytext=(4.5,0.32), textcoords='data',fontsize=lfs,
            arrowprops=dict(facecolor='black',arrowstyle='->'),
            horizontalalignment='right', verticalalignment='top',
            )

ax1.annotate('Jupiter Encounter',xy=(5.4,0.09),xycoords='data',
            xytext=(12,0.42), textcoords='data',fontsize=lfs,
            arrowprops=dict(facecolor='black',arrowstyle='->'),
            horizontalalignment='right', verticalalignment='top',
            )

ax1.annotate('Cover Ejected',xy=(3.2,0.03),xycoords='data',
            xytext=(1.8,-0.25), textcoords='data',fontsize=lfs,
            arrowprops=dict(facecolor='black',arrowstyle='->'),
            horizontalalignment='left', verticalalignment='bottom',
            )

leaderlen = 5
ifield = 1

for i in range(0,np.size(nhpos['fieldflag'])):
    if nhpos['fieldflag'][i,0] == 1:# & ifield < 2 :
        print str(nhpos['fieldflag'][i,0])
        print str(nhpos['elat'][i]) + ':' + str(nhpos['elong'][i])

        thisx = 20 - np.squeeze(nhpos['r'][i])
        thisr = np.sqrt(nhpos['elat'][i]**2 + thisx**2)
        thistheta = nhpos['elat'][i]*np.pi/180
        
        ax2.arrow(np.squeeze(nhpos['r'][i]),\
                  np.squeeze(nhpos['z'][i])/2.6,\
                  np.squeeze(leaderlen * np.cos(thistheta)),\
                  np.squeeze(leaderlen * np.sin(thistheta)),\
                  linestyle='-',color='red',head_width=0.5, head_length=0.5)

        #if np.squeeze(nhpos['r'][i])+\
        #   1.0*np.squeeze(leaderlen * \
        #                  np.cos(nhpos['elong'][i]*np.pi/180)) < 0:
        padderx = 1.2
        paddery = 1.15
        #else:
        #    padderx = 1.25
        #    paddery = 1.4
        
        ax2.annotate('Field '+ str(ifield),\
                     xy=(np.squeeze(nhpos['r'][i]),\
                         np.squeeze(nhpos['z'][i]/2.6)),\
                     xycoords='data',\
                     xytext=(np.squeeze(nhpos['r'][i])+\
                             padderx*np.squeeze(leaderlen * \
                                        np.cos(thistheta)),\
                             np.squeeze(nhpos['z'][i]/2.6)+\
                             paddery*np.squeeze(leaderlen * \
                                                np.sin(thistheta))),\
                     textcoords='data',fontsize=lfs,\
                     horizontalalignment='left',verticalalignment='center')
        
        ax2.plot(nhpos['r'][i],nhpos['z'][i]/2.6,color='red',marker='o')
        ifield = ifield + 1



#plt.show()
plt.tight_layout(w_pad=5.5)
pdf = PdfPages('nh_plot_orbits.pdf')
pdf.savefig()
pdf.close()
