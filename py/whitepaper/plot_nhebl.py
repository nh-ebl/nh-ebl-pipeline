##################################################################
##
##  plot_nhebl.py
##  Nov 9 2015
##  Mike Zemcov
##  This function makes a plot of the various components of a
##  background measurement as made by New Horizons.
##
##################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

ax = plt.subplot(111)

xmin = 0.3
xmax = 3.1
ymin = 1e-2
ymax = 3e3

###################################################################
##  ZL

text_file = open('nh_lookup/ZL_Akari.txt')
rows = [[float(x) for x in line.split()[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

akari_l = np.asarray(cols[0])
akari_y = np.asarray(cols[1])
akari_e = np.asarray(cols[2])

#plt.errorbar(akari_l,akari_y,yerr=akari_e,linestyle='',marker='o',\
#             color='darkred')

text_file = open('nh_lookup/ZL_CIBER.txt')
rows = [[float(x) for x in line.split()[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

ciber_l = np.asarray(cols[0])
ciber_y = np.asarray(cols[1])
ciber_e = np.asarray(cols[2])

#plt.errorbar(ciber_l,ciber_y,yerr=ciber_e,linestyle='',marker='s',\
#             color='darkred')

text_file = open('nh_lookup/ZL_DIRBE.txt')
rows = [[float(x) for x in line.split()[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

dirbe_l = np.asarray(cols[0])
dirbe_y = np.asarray(cols[1])
dirbe_e = np.asarray(cols[2])

#plt.errorbar(dirbe_l,dirbe_y,yerr=dirbe_e,linestyle='',marker='v',\
#             color='darkred')

text_file = open('nh_lookup/zodi_extrap.csv')
rows = [[float(x) for x in line.split()[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

zodi_l = np.asarray(cols[0])
zodi_a = np.asarray(cols[1])
zodi_b = np.asarray(cols[2])
zodi_c = np.asarray(cols[3])
zodi_d = np.asarray(cols[4])
zodi_e = np.asarray(cols[5])
zodi_f = np.asarray(cols[6])

plt.plot(zodi_l,zodi_a,linestyle='-',color='grey')

plt.plot(zodi_l,zodi_f,linestyle='-',color='grey')

plt.plot(zodi_l,zodi_f/50,linestyle='-',color='black',linewidth=1.5)

ax.annotate(r'Zodiacal Light at 1 AU',xy=(0.0,0.0),xytext=(0.68,0.75),xycoords='figure fraction',textcoords='figure fraction',rotation=360-17,color='grey',fontsize=10)

ax.annotate(r'Zodiacal Light at 10 AU',xy=(0.0,0.0),xytext=(0.68,0.385),xycoords='figure fraction',textcoords='figure fraction',rotation=360-17,color='grey',fontsize=10)

ax.annotate(r'Min. Zodiacal Light at 45 AU',xy=(0.0,0.0),xytext=(0.5,0.195),xycoords='figure fraction',textcoords='figure fraction',rotation=360-13,color='black',fontsize=10)

###################################################################
##  DGL

ohmb = 2.

#Ienaka et al 2013
ienaka_x = np.array([0.44, 0.49, 0.55, 0.65])
ienaka_y = ohmb * np.array([10.9, 14.071, 21.803, 15.6814])
ienaka_e = np.array([0.7, 0.6, 1.6, 0.92])

plt.errorbar(ienaka_x,ienaka_y,yerr=ienaka_e,linestyle='',marker='o',\
             color='forestgreen')

#Matsuoka et al 2011
matsuoka_x = np.array([0.44, 0.64])
matsuoka_y = ohmb * np.array([14.308, 21.447])
matsuoka_e = np.array([0.68, 0.468])

plt.errorbar(matsuoka_x,matsuoka_y,yerr=matsuoka_e,linestyle='',marker='s',\
             color='forestgreen')

#witt et al 2008
witt_x = np.array([0.46, 0.53, 0.63, 0.83])
witt_y = ohmb * np.array([13.686, 18.666, 19.986, 9.754])
witt_e = np.array([0.647, 2.26, 3.31, 2.525])

plt.errorbar(witt_x,witt_y,yerr=witt_e,linestyle='',marker='v',\
             color='forestgreen')

text_file = open('nh_lookup/dgl_spec.txt')
rows = [[float(x) for x in line.split()[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

cdgl_l = np.asarray(cols[0])/1000.
cdgl_y = ohmb * np.asarray(cols[1])
cdgl_e = np.asarray(cols[2])

plt.errorbar(cdgl_l,cdgl_y,yerr=cdgl_e,linestyle='',marker='o',color='forestgreen')

text_file = open('nh_lookup/zda04_bc03.dat')
rows = [[float(x) for x in line.split()[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

brandt_l = np.asarray(cols[0]) / 10000.
brandt_y = ohmb * 1.73*np.asarray(cols[1]) 

plt.plot(brandt_l,brandt_y,color='forestgreen')

ax.annotate(r'High Latitude Diffuse Galactic Light',xy=(0.0,0.0),xytext=(0.24,0.65),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='forestgreen',fontsize=10)

###################################################################
##  Gamma Ray Background

text_file = open('nh_lookup/aharonian2014_up.txt')
rows = [[float(x) for x in line.split(',')[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

au_l = np.asarray(cols[0])
au_y = np.asarray(cols[1])

text_file = open('nh_lookup/aharonian2014_dn.txt')
rows = [[float(x) for x in line.split(',')[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

ad_l = np.asarray(cols[0])
ad_y = np.asarray(cols[1])

a_l = np.hstack((au_l,ad_l[::-1]))
a_y = np.hstack((au_y,ad_y[::-1]))

ax.fill_between(a_l,a_y,color='gold',alpha=1)




###################################################################
## IGL

lambda_ebl = np.array([0.25,0.35,0.45,0.545,0.641,0.798,1.22,1.63,2.19,3.6,4.5])
mlimall = np.array([1.9,3.0,4.5,5.65,6.65,7.97,9.60,9.34,8.09,4.87,3.28])
mlimallup = np.array([0.4,0.7,1.5,1.73,1.82,2.01,2.40,2.59,2.52,1.72,1.21])
mlimalldn = np.array([0.4,0.5,0.6,0.85,0.92,1.06,1.28,1.29,1.14,0.71,0.49])

e_l = np.hstack((lambda_ebl,lambda_ebl[::-1]))
e_yu = mlimallup+mlimall
e_yd = mlimall-mlimalldn
e_y = np.hstack((e_yu,e_yd[::-1]))

ax.fill_between(e_l,e_y,color='blue',alpha=0.3)

####
## Madau & Pozzetti 2000 (HST UBVIJHK counts)
## http://adsabs.harvard.edu/abs/2000MNRAS.312L...9M
lmadau=[0.36,0.45,0.67,0.81,1.1,1.6,2.2]
fmadau=[2.87,4.57,6.74,8.04,9.71,9.02,7.92]
errmadauup=[0.58,0.73,1.25,1.62,3.00,2.62,2.04]
errmadaudn=[0.42,0.47,0.94,0.92,1.90,1.68,1.21]

errmadau = np.vstack((errmadaudn,errmadauup))

plt.errorbar(lmadau,fmadau,yerr=errmadau,linestyle='',marker='o',color='blue')

####
## Totani (2001) - Subaru
## http://adsabs.harvard.edu/abs/2001ApJ...550L.137T
ltot=[0.3,0.45,0.61,0.81,1.25,2.2]
ftot=[2.7,4.4,6.0,8.1,10.9,8.3]
errtot=[0.3,0.4,0.6,0.8,1.1,0.8]

plt.errorbar(ltot,ftot,yerr=errtot,linestyle='',marker='o',color='blue')

####
## Fazio (2004) - Spitzer 3.5 um, best expressed as an upper limt
## http://adsabs.harvard.edu/abs/2004ApJS..154...39F
lfaz=[3.6]
fazkjy=[5.4]
ferr = 2.0
fup = np.array([1,0]*5)

plt.errorbar(lfaz,fazkjy,yerr=ferr,uplims=True,linestyle='',marker='o',\
             color='blue')

ax.annotate(r'IGL',xy=(0.0,0.0),xytext=(0.325,0.523),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='blue',fontsize=10)

ax.annotate(r'HESS $\gamma$-ray EBL Estimate',xy=(0.0,0.0),xytext=(0.67,0.531),xycoords='figure fraction',textcoords='figure fraction',rotation=360-5,color='gold',fontsize=10)







####################################################################
##  NH LORRI


ax.errorbar(np.sqrt(0.350*0.850),1.4,yerr=0.6,uplims=True,marker='.',markersize=0,linestyle='',color='red',linewidth=2)

plt.plot(np.asarray([0.350,0.850]),np.asarray([1.4,1.4]),linestyle='-',linewidth=3,color='red')

ax.annotate(r'LORRI Sensitivity in 1 hr',xy=(0.0,0.0),xytext=(0.22,0.42),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='red',fontsize=10)

ax.errorbar(np.sqrt(0.350*0.850),0.1*0.8,yerr=0.04,uplims=True,marker='.',markersize=0,linestyle='',color='deeppink',linewidth=1)

plt.plot(np.asarray([0.350,0.850]),0.1*np.asarray([0.8,0.8]),linestyle='-',linewidth=2,color='deeppink')

ax.annotate(r'LORRI Systematic Limit',xy=(0.0,0.0),xytext=(0.22,0.24),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='deeppink',fontsize=10)


####################################################################
##  NH LEISA

#ax.errorbar(np.sqrt(0.350*0.850),1.4,yerr=0.6,uplims=True,marker='.',markersize=0,linestyle='',color='red',linewidth=2)

leisa_lambda = np.logspace(np.log10(1.25+0.125/2),np.log10(2.5-0.25/2),num=10)

leisa_sens = np.linspace(,num=10)

print leisa_lambda

#plt.plot(leisa_lambda,

'''
plt.plot(np.asarray([0.350,0.850]),np.asarray([1.4,1.4]),linestyle='-',linewidth=3,color='red')

ax.annotate(r'LORRI Sensitivity in 1 hr',xy=(0.0,0.0),xytext=(0.22,0.42),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='red',fontsize=10)

ax.errorbar(np.sqrt(0.350*0.850),0.1*0.8,yerr=0.04,uplims=True,marker='.',markersize=0,linestyle='',color='deeppink',linewidth=1)

plt.plot(np.asarray([0.350,0.850]),0.1*np.asarray([0.8,0.8]),linestyle='-',linewidth=2,color='deeppink')

ax.annotate(r'LORRI Systematic Limit',xy=(0.0,0.0),xytext=(0.22,0.24),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='deeppink',fontsize=10)
'''
'''
cibert_darkcurrent = 0.05 # e-/s

cibert_photocurrent = np.array([9.5,6.8,8.1,7.8,7.7,3.8]) #e-/s
cibert_photonoise = np.array([0.43,0.35,0.41,0.40,0.41,0.35])

cibert_calfac = cibert_s / cibert_photonoise

print cibert_calfac

cibert_darkeff = cibert_calfac * cibert_darkcurrent

print cibert_darkeff

for iis in range(len(cibert_l)):
    ciber_lp = np.array([cibert_l[iis] - cibert_dlh[iis],\
                         cibert_l[iis] + cibert_dlh[iis]])
    ciber_spp = np.array([cibert_darkeff[iis],cibert_darkeff[iis]])
    ax.errorbar(cibert_l[iis],ciber_spp[0],yerr=0.5*ciber_spp[0],uplims=True,marker='.',markersize=0,linestyle='',color='deeppink')
    if iis == 0:
        ax.plot(ciber_lp,ciber_spp,color='deeppink',label='CIBER-2 Dark Current, 100% Error',linewidth=2)
    else:
        ax.plot(ciber_lp,ciber_spp,color='deeppink',linewidth=2)
    if iis == 0:
        ax.plot(ciber_lp,0.01*ciber_spp,color='deeppink',label='CIBER-2 Dark Current, 1% Error',linewidth=2,linestyle='--')
    else:
        ax.plot(ciber_lp,0.01*ciber_spp,color='deeppink',linewidth=2,linestyle='--')
'''
        
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$\lambda$ ($\mu$m)')
ax.set_ylabel(r'$\lambda I_{\lambda}$ (nW m$^{-2}$ sr$^{-1}$)')
ax.set_yscale('log')
ax.set_xscale('log')

plt.legend(fontsize=7,loc=1,handlelength=4)

x=[0.5,0.7,1.0,1.5,2.0,3.0,4.0,5.0]
labels=['0.5','0.7','1.0','1.5','2.0','3.0','4.0','5.0']

#plt.xticks(x, labels)

#plt.show()
plt.savefig('plot_nhebl.png',dpi=200)

plt.close()


   









# return
