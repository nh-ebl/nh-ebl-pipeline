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

fig=plt.figure(figsize=(6.5,4.5))

ax = plt.subplot(111)

xmin = 0.38
xmax = 2.8
ymin = 1e-1
ymax = 1e3

if 0:
    ###################################################################
    ##  ZL

    text_file = open('lookup/ZL_Akari.txt')
    rows = [[float(x) for x in line.split()[:]] for line in text_file]
    cols = [list(col) for col in zip(*rows)]
    text_file.close

    akari_l = np.asarray(cols[0])
    akari_y = np.asarray(cols[1])
    akari_e = np.asarray(cols[2])

    #plt.errorbar(akari_l,akari_y,yerr=akari_e,linestyle='',marker='o',\
        #             color='darkred')

    text_file = open('lookup/ZL_CIBER.txt')
    rows = [[float(x) for x in line.split()[:]] for line in text_file]
    cols = [list(col) for col in zip(*rows)]
    text_file.close

    ciber_l = np.asarray(cols[0])
    ciber_y = np.asarray(cols[1])
    ciber_e = np.asarray(cols[2])

    #plt.errorbar(ciber_l,ciber_y,yerr=ciber_e,linestyle='',marker='s',\
        #             color='darkred')

    text_file = open('lookup/ZL_DIRBE.txt')
    rows = [[float(x) for x in line.split()[:]] for line in text_file]
    cols = [list(col) for col in zip(*rows)]
    text_file.close

    dirbe_l = np.asarray(cols[0])
    dirbe_y = np.asarray(cols[1])
    dirbe_e = np.asarray(cols[2])

    #plt.errorbar(dirbe_l,dirbe_y,yerr=dirbe_e,linestyle='',marker='v',\
        #             color='darkred')

text_file = open('lookup/zodi_extrap.csv')
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

plt.plot(zodi_l,0.2*zodi_f,linestyle='-',color='grey')

#plt.plot(zodi_l,zodi_f/50,linestyle='-',color='black',linewidth=1.5)

ax.annotate(r'IPD Light at 1 AU',xy=(0.0,0.0),xytext=(0.6,0.87),xycoords='figure fraction',textcoords='figure fraction',rotation=360-15,color='grey',fontsize=10)
    
ax.annotate(r'IPD Light at 10 AU',xy=(0.0,0.0),xytext=(0.6,0.24),xycoords='figure fraction',textcoords='figure fraction',rotation=360-14,color='grey',fontsize=10)

#ax.annotate(r'Min. Zodiacal Light at 45 AU',xy=(0.0,0.0),xytext=(0.5,0.195),xycoords='figure fraction',textcoords='figure fraction',rotation=360-13,color='black',fontsize=10)

if 0:
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

    plt.errorbar(matsuoka_x,matsuoka_y,yerr=matsuoka_e,linestyle='',\
                 marker='s',color='forestgreen')

    #witt et al 2008
    witt_x = np.array([0.46, 0.53, 0.63, 0.83])
    witt_y = ohmb * np.array([13.686, 18.666, 19.986, 9.754])
    witt_e = np.array([0.647, 2.26, 3.31, 2.525])

    plt.errorbar(witt_x,witt_y,yerr=witt_e,linestyle='',marker='v',\
                 color='forestgreen')

    text_file = open('lookup/dgl_spec.txt')
    rows = [[float(x) for x in line.split()[:]] for line in text_file]
    cols = [list(col) for col in zip(*rows)]
    text_file.close

    cdgl_l = np.asarray(cols[0])/1000.
    cdgl_y = ohmb * np.asarray(cols[1])
    cdgl_e = np.asarray(cols[2])

    plt.errorbar(cdgl_l,cdgl_y,yerr=cdgl_e,linestyle='',marker='o',color='forestgreen')

    text_file = open('lookup/zda04_bc03.dat')
    rows = [[float(x) for x in line.split()[:]] for line in text_file]
    cols = [list(col) for col in zip(*rows)]
    text_file.close

    brandt_l = np.asarray(cols[0]) / 10000.
    brandt_y = ohmb * 1.73*np.asarray(cols[1]) 

    plt.plot(brandt_l,brandt_y,color='forestgreen')

    ax.annotate(r'High Latitude Diffuse Galactic Light',xy=(0.0,0.0),xytext=(0.24,0.65),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='forestgreen',fontsize=10)
    

###################################################################
##  Gamma Ray Background

text_file = open('lookup/aharonian2014_up.txt')
rows = [[float(x) for x in line.split(',')[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

au_l = np.asarray(cols[0])
au_y = np.asarray(cols[1])

text_file = open('lookup/aharonian2014_dn.txt')
rows = [[float(x) for x in line.split(',')[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

ad_l = np.asarray(cols[0])
ad_y = np.asarray(cols[1])

a_l = np.hstack((ad_l[::-1],au_l))
a_y = np.hstack((ad_y[::-1],au_y))
#a_y[a_y.size-1] = 3.0

ax.fill_between(a_l,a_y,color='silver',alpha=1,interpolate=True)




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

#ax.fill_between(e_l,e_y,color='blue',alpha=0.3)

####
## Madau & Pozzetti 2000 (HST UBVIJHK counts)
## http://adsabs.harvard.edu/abs/2000MNRAS.312L...9M
lmadau=0.99*np.asarray([0.36,0.45,0.67,0.81,1.1,1.6,2.2])
fmadau=[2.87,4.57,6.74,8.04,9.71,9.02,7.92]
errmadauup=[0.58,0.73,1.25,1.62,3.00,2.62,2.04]
errmadaudn=[0.42,0.47,0.94,0.92,1.90,1.68,1.21]

errmadau = np.vstack((errmadaudn,errmadauup))

plt.errorbar(lmadau,fmadau,yerr=errmadau,linestyle='',marker='^',markerfacecolor='none',color='black')
#ax.annotate(r'HDF/SDF',xy=(0.0,0.0),xytext=(0.365,0.455),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='black',fontsize=10)

####
## Totani (2001) - Subaru
## http://adsabs.harvard.edu/abs/2001ApJ...550L.137T
ltot=1.01*np.asarray([0.3,0.45,0.61,0.81,1.25,2.2])
ftot=[2.7,4.4,6.0,8.1,10.9,8.3]
errtot=[0.3,0.4,0.6,0.8,1.1,0.8]

plt.errorbar(ltot,ftot,yerr=errtot,linestyle='',marker='v',color='black',markerfacecolor='none')


####
## Mattila et al 2003 (HST)

#lmattila = np.asarray([0.3,0.555,0.814])
#errmattilaup = np.asarray([9.1,20.7,14.3]) * lmattila * 10
#errmattiladn = np.asarray([7.4,15.4,10.6])*lmattila*10
#fmattila = (errmattilaup + errmattiladn)/2.
#errmattilaup = errmattilaup - fmattila
#errmattiladn = fmattila - errmattiladn
#ax.errorbar(lmattila,fmattila,yerr=[errmattiladn,errmattilaup],marker='^',markerfacecolor='none',color='darkgrey',linestyle='none')


####
## Fazio (2004) - Spitzer 3.5 um, best expressed as an upper limt
## http://adsabs.harvard.edu/abs/2004ApJS..154...39F
lfaz=[3.6]
fazkjy=[5.4]
ferr = 2.0
fup = np.array([1,0]*5)

plt.errorbar(lfaz,fazkjy,yerr=ferr,uplims=True,linestyle='',marker='o',\
             color='blue')

# Gardner 2000 - http://adsabs.harvard.edu/abs/2000ApJ...542L..79G
ax.errorbar([0.2365],[3.6],yerr=[[0.5],[0.7]],marker='s',color='black',markerfacecolor='none',linestyle='')
#ax.annotate(r'STIS',xy=(0.0,0.0),xytext=(0.165,0.4),xycoords='figure fraction',textcoords='figure fraction',rotation=360+0,color='black',fontsize=10)

# Keenan 2010
ax.errorbar([0.99*1.25],[11.7],yerr=[[2.6],[5.6]],marker='>',color='black',markerfacecolor='none',linestyle='')
ax.errorbar([0.99*1.6],[11.5],yerr=[[1.5],[4.5]],marker='>',color='black',markerfacecolor='none',linestyle='')

# Matsuura 2017
mat_wl=np.asarray([1.05,1.11,1.18,1.25,1.33,1.42,1.51,1.60,1.70])	
mat_il=np.asarray([7.91,15.39,12.86,20.13,24.27,28.68,23.79,27.05,24.82])	
mat_er=np.asarray([3.83,3.44,3.18,3.05,2.84,2.99,3.52,3.79,3.016])

ax.errorbar(mat_wl,mat_il,yerr=mat_er,color='grey',marker='p',linestyle='')
    
# dark cloud Mattila 2017
ax.errorbar([0.400],[11.6],yerr=[[4.4],[4.4]],color='darkgrey',linewidth=1,marker='s')
ax.plot([0.430],[20.0],color='darkgrey',marker='s')
ax.errorbar([0.430],[20.0],yerr=[10],uplims=True,color='darkgrey',linewidth=1)
ax.plot([0.520],[23.4],color='darkgrey',marker='s')
ax.errorbar([0.520],[23.4],yerr=[10],uplims=True,color='darkgrey',linewidth=1)



####################################################################

#ax.annotate(r'IGL',xy=(0.0,0.0),xytext=(0.325,0.523),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='blue',fontsize=10)

ax.annotate(r'HESS $\gamma$-rays',xy=(0.0,0.0),xytext=(0.38,0.555),xycoords='figure fraction',textcoords='figure fraction',rotation=360+8,color='black',fontsize=10)


# Bernstein 2007
#ax.errorbar(np.asarray([0.300,0.300]),3*np.asarray([6.0,6.0]),yerr=3*np.asarray([4.0,4.0]),marker='s',color='black',markerfacecolor='none')
#ax.errorbar(np.asarray([0.555,0.555]),5.55*np.asarray([10,10]),yerr=5.55*np.asarray([5.0,5.0]),marker='s',color='black',markerfacecolor='none')
#ax.errorbar(np.asarray([0.814,0.814]),8.14*np.asarray([7.0,7.0]),yerr=8.14*np.asarray([4.0,4.0]),marker='s',color='limegreen')
#ax.annotate(r'WFPC2',xy=(0.0,0.0),xytext=(0.475,0.675),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='black',fontsize=10)

# Cambresy 2001
ax.errorbar(0.97*np.asarray([1.25,1.25]),np.asarray([54.0,54.0]),yerr=np.asarray([16.8,16.8]),color='lightgrey',marker='D')

# Levenson 2007
ax.errorbar(0.99*np.asarray([1.25,1.25]),2.99*np.asarray([8.9,8.9])/1.25,yerr=2.99*np.asarray([6.3,6.3])/1.25,color='lightgrey',marker='D')

# Wright 2001 - http://iopscience.iop.org/article/10.1086/320942/pdf
ax.errorbar(1.01*np.asarray([1.25,1.25]),np.asarray([28.9,28.9]),yerr=np.asarray([16.3,16.3]),color='lightgrey',marker='D')

# Wright 2004 - http://adsabs.harvard.edu/abs/2004NewAR..48..465W
ax.errorbar(1.03*np.asarray([1.25,1.25]),np.asarray([28.0,28.0]),yerr=np.asarray([15,15]),color='lightgrey',marker='D')
ax.annotate(r'DIRBE',xy=(0.0,0.0),xytext=(0.55,0.63),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='black',fontsize=10)


# Matsumoto (2005) - IRTS
# http://adsabs.harvard.edu/abs/2005ApJ...626...31M
mat05nw=[15.5,15.3,13.5,13.1,15.5,14.4,14.6,13.6,16.4,19.7,\
         18.6,19.5,23.7,22.7,24.4,29.7,35.4,39.2,43.2,51.0,\
         58.7,65.9,71.3,70.1]
emat05nw=[3.9,3.6,3.3,3.7,3.0,3.0,3.0,3.0,2.9,3.0,3.2,3.5,4.2,\
          4.5,4.8,5.3,6.0,7.0,7.9,8.8,10.3,11.9,12.8,13.2]
lmat05=[3.98,3.88,3.78,3.68,3.58,3.48,3.38,3.28,3.17,3.07,2.98,\
        2.88,2.54,2.44,2.34,2.24,2.14,2.03,1.93,1.83,1.73,1.63,\
        1.53,1.43]
ax.errorbar(lmat05,mat05nw,yerr=emat05nw,marker='+',color='black',linestyle='')
ax.annotate(r'IRTS',xy=(0.0,0.0),xytext=(0.83,0.67),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='black',fontsize=10)

# Zemcov 2014 - http://adsabs.harvard.edu/abs/2014Sci...346..732Z
ax.errorbar(np.asarray([1.1]),np.asarray([16.7]),yerr=[[4],[5]],color='darkgrey',marker='p')
ax.errorbar(np.asarray([1.6]),np.asarray([20.4]),yerr=[[5.1],[6]],color='darkgrey',marker='p')
#ax.annotate(r'CIBER',xy=(0.0,0.0),xytext=(0.705,0.615),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='black',fontsize=10)


#if 0:
####################################################################
##  NH LORRI

lorri_tint = 3600./30.
ax.errorbar(np.sqrt([0.440*0.870]),10./np.sqrt(lorri_tint),yerr=0.4,uplims=True,marker='.',markersize=0,linestyle='',color='blue',linewidth=1)
    
plt.plot(np.asarray([0.440,0.870]),np.asarray([10,10])/np.sqrt(lorri_tint),linestyle='-',linewidth=1,color='blue')

#ax.annotate(r'Hypothetical LORRI Survey',xy=(0.0,0.0),xytext=(0.36,0.36),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='red',fontsize=8)

if 0:
    ax.errorbar(np.sqrt(0.350*0.850),0.1*0.8,yerr=0.04,uplims=True,marker='.',markersize=0,linestyle='',color='deeppink',linewidth=1)

    plt.plot(np.asarray([0.350,0.850]),0.1*np.asarray([0.8,0.8]),linestyle='-',linewidth=2,color='deeppink')

    ax.annotate(r'LORRI Systematic Limit',xy=(0.0,0.0),xytext=(0.22,0.24),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,color='deeppink',fontsize=10)

# matsuoka et al 2011
ax.errorbar(np.asarray([0.44,0.44]),np.asarray([7.9,7.9]),yerr=np.asarray([4.0,4.0]),marker='o',color='dimgray')
ax.errorbar(np.asarray([0.64,0.64]),np.asarray([7.7,7.7]),yerr=np.asarray([5.8,5.8]),marker='o',color='dimgray')
ax.plot([0.3950,0.4850],[7.9,7.9],color='dimgray',linestyle='-')
ax.plot([0.59,0.69],[7.7,7.7],color='dimgray',linestyle='-')

# Toller 1983
ax.plot([0.44],[19.8],marker='o',color='dimgrey')
ax.errorbar([0.44],[19.8],yerr=6, uplims=True,color='dimgrey')
ax.plot([0.44],[19.8],marker='o',color='dimgrey')

ax.annotate(r'Pioneer',xy=(0.0,0.0),xytext=(0.232,0.49),xycoords='figure fraction',textcoords='figure fraction',rotation=0,color='black',fontsize=10)

# LORRI
#ax.plot(np.asarray([0.440,0.870]),np.asarray([9,9]),color='red',linestyle='-',linewidth=1.)
#eb1=ax.errorbar(np.asarray([0.655,0.655]),np.asarray([9,9]),yerr=np.asarray([7.3,7.3]/np.sqrt(1165/260)),marker='o',color='red',linewidth=1.,markersize=10,markerfacecolor='red',markeredgecolor='none')
#eb1[-1][0].set_linestyle('--')
ax.plot(np.asarray([0.440,0.870]),np.asarray([19.3,19.3]),color='blue',linestyle='--',linewidth=1)
ax.errorbar([np.sqrt(0.44*0.870)],[19.3],yerr=9, uplims=True,color='blue',linewidth=1)
ax.annotate(r'Current LORRI Limit',xy=(0.0,0.0),xytext=(0.28,0.605),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,fontsize=6,color='blue')
#ax.annotate(r'Proposed Study',xy=(0.505,0.525),xytext=(0.57,0.405),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,fontsize=10,color='red',arrowprops=dict(facecolor='red', shrink=0.1))
#ax.plot(np.asarray([0.440,0.870]),np.asarray([19.3,19.3]),color='red',linestyle='-',linewidth=1)
leisa_wl = np.logspace(np.log10(1.25),np.log10(2.5),num=10)
leisa_sens = np.repeat(6e4,10)/np.sqrt(240./10.)/np.sqrt(256.)/np.sqrt(256./240.)
ax.errorbar(np.sqrt(1.25*2.5),leisa_sens[0]/np.sqrt(24.*3600./4.),yerr=2, uplims=True,color='red',linestyle='-')
ax.plot(leisa_wl,leisa_sens/np.sqrt(24.*3600./4.),color='red',linestyle='-',marker='o')
ax.annotate(r'LEISA (1 day, $R=10$)',xy=(0.0,0.0),xytext=(0.64,0.41),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,fontsize=12,color='red')
#ax.plot(np.asarray([0.440,0.870]),np.asarray([19.3,19.3]/np.sqrt(7200./260)),color='red',linestyle='-',linewidth=2)
mvic_tint = 1.*3600.*4.
mvic_wl1 = np.array([0.4,0.975])
mvic_sens1 = np.repeat(7.5e4,2)/np.sqrt(2.*5000.*32.)
ax.errorbar(np.sqrt(mvic_wl1[0]*mvic_wl1[1]),mvic_sens1[0]/np.sqrt(mvic_tint),yerr=0.5, uplims=True,color='green',linestyle='-')
ax.plot(mvic_wl1,mvic_sens1/np.sqrt(mvic_tint),color='green',linestyle='-')
mvic_wl2 = np.array([0.4,0.55])
mvic_sens2 = np.repeat(7.5e4,2)/np.sqrt(5000.*32.)
ax.errorbar(np.sqrt(mvic_wl2[0]*mvic_wl2[1]),mvic_sens2[0]/np.sqrt(mvic_tint),yerr=0.5, uplims=True,color='turquoise',linestyle='-')
ax.plot(mvic_wl2,mvic_sens2/np.sqrt(mvic_tint),color='turquoise',linestyle='-')
mvic_wl3 = np.array([0.540,0.700])
mvic_sens3 = np.repeat(7.5e4,2)/np.sqrt(5000.*32.)
ax.errorbar(np.sqrt(mvic_wl3[0]*mvic_wl3[1]),0.95*mvic_sens3[0]/np.sqrt(mvic_tint),yerr=0.5, uplims=True,color='springgreen',linestyle='-')
ax.plot(mvic_wl3,0.95*mvic_sens3/np.sqrt(mvic_tint),color='springgreen',linestyle='-')
mvic_wl4 = np.array([0.780,0.975])
mvic_sens4 = np.repeat(7.5e4,2)/np.sqrt(5000.*32.)
ax.errorbar(np.sqrt(mvic_wl4[0]*mvic_wl4[1]),mvic_sens4[0]/np.sqrt(mvic_tint),yerr=0.5, uplims=True,color='limegreen',linestyle='-')
ax.plot(mvic_wl4,mvic_sens4/np.sqrt(mvic_tint),color='limegreen',linestyle='-')
mvic_wl5 = np.array([0.860,0.910])
mvic_sens5 = np.repeat(7.5e4,2)/np.sqrt(5000.*32.)
ax.errorbar(np.sqrt(mvic_wl5[0]*mvic_wl5[1]),2*mvic_sens5[0]/np.sqrt(mvic_tint),yerr=0.5, uplims=True,color='yellowgreen',linestyle='-')
ax.plot(mvic_wl5,2*mvic_sens5/np.sqrt(mvic_tint),color='yellowgreen',linestyle='-')
ax.annotate(r'MVIC (1 hr)',xy=(0.0,0.0),xytext=(0.28,0.4),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,fontsize=12,color='green')
ax.annotate(r'LORRI (1 hr)',xy=(0.0,0.0),xytext=(0.345,0.3),xycoords='figure fraction',textcoords='figure fraction',rotation=360-0,fontsize=12,color='blue')


ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.set_xlabel(r'$\lambda$ ($\mu$m)')
ax.set_ylabel(r'$\lambda I_{\lambda}^{\rm COB}$ (nW m$^{-2}$ sr$^{-1}$)')
ax.set_yscale('log')
ax.set_xscale('log')

#plt.legend(fontsize=7,loc=1,handlelength=4)

x=[0.4,0.5,0.6,0.7,0.8,0.9,\
   1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,\
   2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8]
labels=['0.4','','','','','',\
        '1.0','','','','','','','','','',\
        '2.0','','','','','','','','2.8']

plt.xticks(x, labels)

plt.tight_layout()
#plt.show()
plt.savefig('paper_plot_cob.pdf')

plt.close()

### to add - Old Pioneer
   







# return
