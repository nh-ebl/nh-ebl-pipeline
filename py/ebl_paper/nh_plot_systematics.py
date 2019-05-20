import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

fig=plt.figure(figsize=(6.5,5))

ax = plt.subplot(111)

xmin = 0.5
xmax = 7.5
ymin=-7
ymax=15

num=7
width = 0.35
ind = np.arange(num) - width / 2 + 1
err_lo = np.array([0.9,3.8,1.9,2.6,0.1,0.7,9.1])
err_hi = np.array([0.0,3.8,1.9,0.0,0.1,0.7,8.6])
meann = 4.7
mean = meann * np.ones(num)

colors=['grey','blue','cyan','seagreen','green','lawngreen','gold','orange']

ax.plot([1.5,1.5],[ymin,ymax],linestyle='--',color='black')
ax.plot([3.5,3.5],[ymin,ymax],linestyle='--',color='black')

ax.bar(ind,err_lo,width,bottom=mean - err_lo,color=colors)
ax.bar(ind,err_hi,width,bottom=mean,color=colors)
ax.plot([xmin,xmax],[meann,meann],linestyle='-',color='red')
ax.plot([xmin,xmax],[0,0],linestyle=':',color='black')
ax.plot([xmin,xmax],[mean + 7.3,mean + 7.3],linestyle='--',color='red')
ax.plot([xmin,xmax],[mean - 7.3,mean - 7.3],linestyle='--',color='red')

ax.annotate('Inst.',xy=(0.75,-6),xytext=(0.75,-6))
ax.annotate('Calibration',xy=(1.91,-6),xytext=(1.91,-6))
ax.annotate('Astrophysical',xy=(4.7,-6),xytext=(4.7,-6))

ax.set_xticklabels(('', 'Dark Current', 'Photometric Calibration',r'$\Omega_{\rm beam}$ Uncertainty', 'IPD Light','USNO-B1 Uncertainty','Star Sample Variance','DGL Uncertainty',''),rotation='vertical')

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_ylabel(r'$\delta \lambda I_{\lambda}^{\rm COB}$ (nW m$^{-2}$ sr$^{-1}$)')

plt.tight_layout()
#plt.show()
plt.savefig('nh_plot_systematics.pdf')

plt.close()
