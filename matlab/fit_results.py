"""
Plot NH LORRI results.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import h5py


FLG_useError = True #turn off x and y errors being used (just basic linear fit if off)

FLG_debug = False #plots various extra lines wrt to the post-extinction green line


#############################  Preamble  ##################################
plt.rcParams.update({'font.size': 15})

## set up plotting
plt.ioff() #disables the plot from showing, which is fine since it's saved directly anyway
fig=plt.figure(figsize=(6,5))
ax = fig.add_subplot(1,1,1)

## read in Teresa's data file (I manually excised the header line)
## and make variables we will use
text_file = open('fit_info.txt')
rows = [[float(x) for x in line.split('\t')[:]] for line in text_file]
cols = [list(col) for col in zip(*rows)]
text_file.close

lil_opt = np.asarray(cols[0])
lil_100m = np.asarray(cols[1])
# lil_db = np.asarray(cols[2]) #removed since multiplication happens at input into function
lil_sigy = np.asarray(cols[2])
lil_sigx = np.asarray(cols[3])
lil_ext = np.asarray(cols[4])
run_num = np.asarray(cols[5])

# Set y-errors to zero
# if np.all(run_num == 2):
#     lil_sigy = np.zeros(lil_sigy.shape,dtype=lil_sigy.dtype)

# lil_x = lil_100m * lil_db
lil_x = lil_100m #alias, already includes db (dl) now

#Read in plot text file
text_file = open('fit_info_txt.txt')
rows = [[x for x in line.split('\t')[:]] for line in text_file][0]
# cols = [list(col) for col in zip(*rows)]
text_file.close

title = rows[0] #title of data

#Prep a text file to save the printed text
if np.all(run_num == 1):
    textOut_file = open('nh_text_results_2022_'+title.replace('/','')+'.txt','w')
elif np.all(run_num == 2):
    textOut_file = open('nh_text_results_nodb_2022_'+title.replace('/','')+'.txt','w')

## this is a guess at b(lambda) used in the error estimation for the
## zeroth-cut fit.  Results shouldn't be sensitive to it if it's about
## in the right place.
if( title == 'NHI' ):
    slope_guess = (lil_opt[np.max(lil_x) == lil_x][0]-lil_opt[np.min(lil_x) == lil_x][0])/(np.max(lil_x)-np.min(lil_x))
else:
    slope_guess = 7.

## plot the "as given" data in blue
ax.errorbar(lil_x,lil_opt,yerr=lil_sigy,xerr=lil_sigx,linestyle='',marker='o',markerfacecolor='xkcd:dark lilac',markeredgecolor='xkcd:dark lilac',ecolor='xkcd:dark lilac')

#############################  Part 1  ###################################

## do a fit with the zeroth-cut design matrix.  This fit isn't what we
## would report, but gives a good "gut instinct" about how big the COB will be

## first, build the design matrix and noise estimate (see Numerical
## Recipes Chapter 15)

A = np.ones([np.size(lil_opt),2])
Ninv = np.zeros([np.size(lil_opt)])
for ind in range(0,np.size(lil_opt)):
    A[ind,1] = lil_x[ind]
    A[ind,0] = 1.
    Ninv[ind] = 1./(lil_sigy[ind]**2 + slope_guess**2 * lil_sigx[ind]**2)

## now do chisq minimization using linear algebra
At = np.transpose(A)
AtAinv = np.linalg.inv((np.matmul(At*Ninv,A)))
AtNd = np.matmul(At*Ninv,lil_opt)
## m will hold the best fit parameters
if( FLG_useError ):
    m = np.matmul(AtAinv,AtNd)
else:
    m = np.flip(np.polyfit(lil_x, lil_opt, 1))

## now tell me what you got.
if np.all(run_num == 1):
    print('WITH D(B) INCLUDED')
elif np.all(run_num == 2):
    print('WITHOUT D(B) INCLUDED')

print('Zeroth-order cut - informational only. '+title)
print('Zeroth-order cut - informational only. '+title,file=textOut_file)
if( title == 'NHI' ):
    print('Best-fitting slope is: ',m[1],' \pm ',np.sqrt(AtAinv[1,1]),'.')
    print('Best-fitting offset is: ',m[0],' \pm ',np.sqrt(AtAinv[0,0]),'.')
    print('Best-fitting slope is: ',m[1],' \pm ',np.sqrt(AtAinv[1,1]),'.',file=textOut_file)
    print('Best-fitting offset is: ',m[0],' \pm ',np.sqrt(AtAinv[0,0]),'.',file=textOut_file)
else:
    print('Best-fitting slope is: {0:2.2f} \pm {1:2.2f}.'.format(m[1],np.sqrt(AtAinv[1,1])))
    print('Best-fitting offset is: {0:2.2f} \pm {1:2.2f}.'.format(m[0],np.sqrt(AtAinv[0,0])))
    print('Best-fitting slope is: {0:2.2f} \pm {1:2.2f}.'.format(m[1],np.sqrt(AtAinv[1,1])),file=textOut_file)
    print('Best-fitting offset is: {0:2.2f} \pm {1:2.2f}.'.format(m[0],np.sqrt(AtAinv[0,0])),file=textOut_file)

## show the output on a plot with a solid blue line
if( title == 'NHI' ):
    mod_x = np.linspace(0,5.5E20)
else:
    mod_x = np.linspace(0,4.5)

mod_y = m[0] + m[1] * mod_x
COBEst_orig, = ax.plot(mod_x,mod_y,color='xkcd:dark lilac')

#############################  Part 2  ###################################

## now do a fit including offset knowledge of extinctions.  This first
## run is just to get the slope correct for the error estimate, and isn't
## the final answer, so I won't print it out.
## build design matrix

A1 = np.ones([np.size(lil_opt),2])
Ninv1 = np.zeros([np.size(lil_opt)])
for ind in range(0,np.size(lil_opt)):
    A1[ind,1] = lil_x[ind]
    A1[ind,0] = lil_ext[ind]
    Ninv1[ind] = 1./(lil_sigy[ind]**2 + m[1]**2 * lil_sigx[ind]**2)

## do linear algebra
At1 = np.transpose(A1)
AtAinv1 = np.linalg.inv((np.matmul(At1*Ninv1,A1)))
AtNd1 = np.matmul(At1*Ninv1,lil_opt)
m1 = np.matmul(AtAinv1,AtNd1)

## now we have a good estimate for m1, we can do the real fit
## buld the design matrix (again, just to be sure it's correct)
for ind in range(0,np.size(lil_opt)):
    A1[ind,1] = lil_x[ind]
    A1[ind,0] = lil_ext[ind]
    Ninv1[ind] = 1./(lil_sigy[ind]**2 + m1[1]**2 * lil_sigx[ind]**2)
## do the linear algebra
AtAinv1 = np.linalg.inv((np.matmul(At1*Ninv1,A1)))
AtNd1 = np.matmul(At1*Ninv1,lil_opt)
## m1 holds the parameter estimates
if( FLG_useError ):
    m1 = np.matmul(AtAinv1,AtNd1)
else:
    m1 = np.flip(np.polyfit(lil_x, lil_opt, 1))

# Calculate RMSE of fit
guess_y = m[1]*lil_x+m[0]
res2 = (guess_y-lil_opt)**2
rss = np.sum(res2*Ninv)
dof = len(lil_opt) - 2
mse = rss/dof
og_rmse = np.sqrt(mse)
if( title == 'NHI' ):
    print('Original-fit RMSE is: ',og_rmse,'.')
    print('Original-fit RMSE is: ',og_rmse,'.',file=textOut_file)
else:
    print('Original-fit RMSE is: {0:2.2f}.'.format(og_rmse))
    print('Original-fit RMSE is: {0:2.2f}.'.format(og_rmse),file=textOut_file)
ax.fill_between(mod_x, mod_y-og_rmse,mod_y+og_rmse, facecolor='xkcd:dark lilac', alpha=0.5,edgecolor='none')

## now tell me what you got.
print('Extinction-corrected COB estimates (USE THESE NUMBERS). '+title)
print('Extinction-corrected COB estimates (USE THESE NUMBERS). '+title,file=textOut_file)
if( title == 'NHI' ):
    print('Best-fitting slope is: ',m1[1],' \pm ',np.sqrt(AtAinv1[1,1]),'.')
    print('Best-fitting offset is: ',m1[0],' \pm ',np.sqrt(AtAinv1[0,0]),'.')
    print('Best-fitting slope is: ',m1[1],' \pm ',np.sqrt(AtAinv1[1,1]),'.',file=textOut_file)
    print('Best-fitting offset is: ',m1[0],' \pm ',np.sqrt(AtAinv1[0,0]),'.',file=textOut_file)
else:
    print('Best-fitting slope is: {0:2.2f} \pm {1:2.2f}.'.format(m1[1],np.sqrt(AtAinv1[1,1])))
    print('Best-fitting offset is: {0:2.2f} \pm {1:2.2f}.'.format(m1[0],np.sqrt(AtAinv1[0,0])))
    print('Best-fitting slope is: {0:2.2f} \pm {1:2.2f}.'.format(m1[1],np.sqrt(AtAinv1[1,1])),file=textOut_file)
    print('Best-fitting offset is: {0:2.2f} \pm {1:2.2f}.'.format(m1[0],np.sqrt(AtAinv1[0,0])),file=textOut_file)

## plot it as a red dashed line
mod_y1 = m1[0] + m1[1] * mod_x
COBEst_adjusted, = ax.plot(mod_x,mod_y1,color='xkcd:dull green')
#############################  Part 3  ###################################

## ok, this next part is a check - it should be the case that if I:
## 1) take the input numbers for (e*COB+DGL).
## 2) Fit for COB.  (These two steps are done above).
## 3) Subtract that estimate from the input, and then add back COB/e
##     in each field.
## 4) Fit that to a line.
## I should have corrected the COB estimate for extinction, at
## least to 1st order.  So this gives a better estimate than the
## first guess, and it would be nice to see it to compare against
## the design-matrix method used above.

## make a new estimate for where the data points should be based on
## our guess at the COB value from step 1
lil_optp = lil_opt - m[0] + m[0] / lil_ext
lil_residual = lil_optp - ( m1[1]*lil_x + m1[0] )
lil_residual_prext = lil_opt - ( m[1]*lil_x + m[0] )
with h5py.File('fit_results_data_'+title.replace('/','')+'.h5', 'w') as fandle:
    fandle.create_dataset('lil_optp', lil_optp.shape, data=lil_optp)
    fandle.create_dataset('lil_residual', lil_residual.shape, data=lil_residual)
    fandle.create_dataset('lil_opt', lil_opt.shape, data=lil_opt)
    fandle.create_dataset('lil_residual_prext', lil_residual_prext.shape, data=lil_residual_prext)
    fandle.create_dataset('lil_x', lil_x.shape, data=lil_x)

## make new error estimate using a better guess for the slope (details
## of which shouldn't really matter here)
for ind in range(0,np.size(lil_opt)):
    Ninv[ind] = 1./(lil_sigy[ind]**2 + m1[1]**2 * lil_sigx[ind]**2)

## now do linear algebra again to get parameter estimates
AtAinvp = np.linalg.inv((np.matmul(At*Ninv,A)))
AtNdp = np.matmul(At*Ninv,lil_optp)
# parameter estimates are in mp
if( FLG_useError ):
    mp = np.matmul(AtAinvp,AtNdp)
else:
    mp = np.flip(np.polyfit(lil_x, lil_optp, 1))

# Calculate RMSE of fit
guess_y = m1[1]*lil_x+m1[0]
res2 = (guess_y-lil_optp)**2
rss = np.sum(res2*Ninv1)
dof = len(lil_opt) - 2
mse = rss/dof
adj_rmse = np.sqrt(mse)
if( title == 'NHI' ):
    print('Adjusted-fit RMSE is: ',adj_rmse,'.')
    print('Adjusted-fit RMSE is: ',adj_rmse,'.',file=textOut_file)
else:
    print('Adjusted-fit RMSE is: {0:2.2f}.'.format(adj_rmse))
    print('Adjusted-fit RMSE is: {0:2.2f}.'.format(adj_rmse),file=textOut_file)

ax.fill_between(mod_x, mod_y1-adj_rmse,mod_y1+adj_rmse, facecolor='xkcd:dull green', alpha=0.5,edgecolor='none')

## now tell me what you got.
print('Extinction-bumped COB estimates (check these against previous estimates). '+title)
print('Extinction-bumped COB estimates (check these against previous estimates). '+title,file=textOut_file)
if( title == 'NHI' ):
    print('Best-fitting slope is: ',mp[1],' \pm ',np.sqrt(AtAinvp[1,1]),'.')
    print('Best-fitting offset is: ',mp[0],' \pm ',np.sqrt(AtAinvp[0,0]),'.')
    print('Best-fitting slope is: ',mp[1],' \pm ',np.sqrt(AtAinvp[1,1]),'.',file=textOut_file)
    print('Best-fitting offset is: ',mp[0],' \pm ',np.sqrt(AtAinvp[0,0]),'.',file=textOut_file)
else:
    print('Best-fitting slope is: {0:2.2f} \pm {1:2.2f}.'.format(mp[1],np.sqrt(AtAinvp[1,1])))
    print('Best-fitting offset is: {0:2.2f} \pm {1:2.2f}.'.format(mp[0],np.sqrt(AtAinvp[0,0])))
    print('Best-fitting slope is: {0:2.2f} \pm {1:2.2f}.'.format(mp[1],np.sqrt(AtAinvp[1,1])),file=textOut_file)
    print('Best-fitting offset is: {0:2.2f} \pm {1:2.2f}.'.format(mp[0],np.sqrt(AtAinvp[0,0])),file=textOut_file)

## subtlety!! - I'd like to plot these corrected points against
## the model derived in part 2 since they should be (roughly) comparable.
## I'm plotting them as open face to highlight that they've been
## monkeyed with ie by trying to correct them for the COB extinction     
ax.errorbar(lil_x,lil_optp,yerr=lil_sigy,xerr=lil_sigx,linestyle='',marker='o',markerfacecolor='white',markeredgecolor='xkcd:dull green',ecolor='xkcd:dull green')

#############################  Postscript  ###################################
if( (title != 'NHI') & (FLG_debug == True) ):
    ax.plot(mod_x, 6*mod_x+18,color='xkcd:cyan')
    ax.plot(mod_x, mp[0] + mp[1] * mod_x, color='xkcd:red')
    mp_direct = np.polyfit(lil_x, lil_optp, 1)
    ax.plot(mod_x, mp_direct[1] + mp_direct[0] * mod_x, color='xkcd:orange')
## plotting stuff - I've commented out stuff that might be useful
## but I didn't need in the first instance.
    
#ax.set_ylim([0,120])
if( title == 'NHI' ):
    ax.set_xlim([0.25E20,5.5E20])
else:
    ax.set_xlim([0,4.5])

if np.all(run_num == 1):
    if( title == 'NHI' ):
        ax.set_xlabel(r'NHI * d(b) [cm$^{-2}$]')
    else:
        ax.set_xlabel(r'$\nu I_{\nu}^{\rm 100 \mu m}$ * d(b) [MJy sr$^{-1}$]')
elif np.all(run_num == 2):
    if( title == 'NHI' ):
        ax.set_xlabel(r'NHI [cm$^{-2}$]')
    else:
        ax.set_xlabel(r'$\nu I_{\nu}^{\rm 100 \mu m}$ [MJy sr$^{-1}$]')
    
ax.set_ylabel(r'$\lambda I_{\lambda}^{\rm EBL+DGL}$ [nW m$^{-2}$ sr$^{-1}$]')
ax.set_title(title)
# ax.legend([COBEst_orig,COBEst_adjusted], ['Original COB Estimate','Extinction-Adjusted COB Estimate'],loc='upper left')
ax.legend([COBEst_adjusted,COBEst_orig], ['Correlative COB','Pre-Extinction Fit'],loc='upper left')

if np.all(run_num == 1):
    plt.tight_layout()
    plt.show(block=False)
    pdf = PdfPages('nh_plot_results_2022_'+title.replace('/','')+'.pdf')
    pdf.savefig(bbox_inches='tight')
    pdf.close()
    plt.close()
elif np.all(run_num == 2):
    plt.tight_layout()
    plt.show(block=False)
    pdf = PdfPages('nh_plot_results_nodb_2022_'+title.replace('/','')+'.pdf')
    pdf.savefig(bbox_inches='tight')
    pdf.close()
    plt.close()
