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

darktemp = loadmat('../../scratch/nh_dark_analysis_fig1.mat')

fig=plt.figure(figsize=(10.9,5))


ax = fig.add_subplot(1,2,1)

p3,=ax.semilogx(np.squeeze(darktemp['mydates']),np.squeeze(darktemp['myfunc']),linestyle='-',color='skyblue',linewidth=5)

p1,=ax.semilogx(darktemp['darkdate'],darktemp['darktemp'],marker='o',linestyle='',color='blue')
p2,=ax.semilogx(darktemp['lightdate'],darktemp['lighttemp'],marker='o',linestyle='',color='red')

ax.semilogx(np.squeeze(darktemp['cover']),[190,225],linestyle='--',color='black')

ax.annotate('Cover Ejected',xy=(220,202),xycoords='data',\
            xytext=(220,202),textcoords='data',\
            horizontalalignment='right',verticalalignment='top',\
            rotation=90)

ax.set_ylim([190,225])
ax.set_xlim([50,3000])
ax.set_xlabel('Time from Launch (days)')
ax.set_ylabel('Detector Temperature (K)')

ax.legend([p1,p2,p3],['Cover-on Data','Cover-off Data','Exponential Decrease'],\
          frameon=False,numpoints=1,fontsize=12)

dref = loadmat('../../scratch/nh_dark_analysis_fig6.mat')

axa = fig.add_subplot(1,2,2)

q3,=axa.plot(np.squeeze(dref['tccd']),np.squeeze(dref['modelone']),linestyle='--',color='midnightblue',linewidth=1)
q4,=axa.plot(np.squeeze(dref['tccd']),np.squeeze(dref['modeltwo']),linestyle='-',color='midnightblue',linewidth=1)

axa.errorbar(np.squeeze(dref['lighttempm']),np.squeeze(dref['lightrefm'][:,0]),yerr=np.squeeze(dref['lightrefm'][:,1]),xerr=0.15*np.ones(4),marker='o',color='red',linestyle='')
axa.errorbar(np.squeeze(dref['darktempm']),np.squeeze(dref['darkrefm'][:,0]),yerr=np.squeeze(dref['darkrefm'][:,1]),xerr=0.15*np.ones(4),marker='o',color='blue',linestyle='')

axa.set_xlim([190,225])
axa.set_ylim([539,549])

axa.plot([190,225],[np.squeeze(dref['meanvref']),np.squeeze(dref['meanvref'])],\
         linestyle='--',color='red')
axa.plot([190,225],[np.squeeze(dref['meanvref']+dref['sigvref']),\
                    np.squeeze(dref['meanvref']+dref['sigvref'])],\
         linestyle=':',color='red')
axa.plot([190,225],[np.squeeze(dref['meanvref']-dref['sigvref']),\
                    np.squeeze(dref['meanvref']-dref['sigvref'])],\
         linestyle=':',color='red')        

axa.set_xlabel('Detector Temperature (K)')
axa.set_ylabel('Mean Value of Reference Pixels (DN)')

axa.legend([p1,p2,q4,q3],['Cover-on Data','Cover-off Data','Model Fit',\
                          'Expected Performance'],\
frameon=False,numpoints=1,fontsize=12,loc=2)


plt.tight_layout(w_pad=5.5)

#plt.show()
pdf = PdfPages('nh_plot_reference.pdf')
pdf.savefig()
pdf.close()
