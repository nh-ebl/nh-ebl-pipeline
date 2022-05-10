# Before running this code you MUST edit /home/symons/anaconda3/lib/python3.7/site-packages/astrobase/services/trilegal.py
# ADD to TRILEGAL_FILTER_SYSTEMS dict:
#     'gaiaDR2': {
#         'desc': "Gaia's DR2 G, G_BP, G_RP (Vegamags, Gaia passbands from Evans et al. 2018)",
#         'table': 'tab_mag_odfnew/tab_mag_gaiaDR2.dat'
#     },
# COMMENT OUT extinction_info and Av_infinity lines w/ try/except ect. in query_galcoords function and replace them with:
#     #instead, use website default
#     Av_infinity = 0.0378
#     inputparams['extinction_infty'] = '%.5f' % Av_infinity

import os
import gzip
import shutil
import astrobase.services.trilegal

dir_base = '/data/symons/nh_data_lauer/trilegal/gaia' #base directory for outputs

# lauer fields
field_DEC = [-17.7780, -41.6360, -50.7529, -0.5183, 0.2915, 36.2331, 35.2979] #DEC values
field_RA = [1.7790, 348.0611, 33.4069, 358.2428, 0.8066, 224.9875, 226.4865] #RA values
field_fieldNums = [i for i in range(1,len(field_DEC)+1)] #field num array, make a list if not going 1 to lenDEC
field_numPerField = 10 #num times to run the same field

field_numFields = len(field_DEC) #get len now
if( (len(field_DEC) != len(field_RA)) | (len(field_DEC) != len(field_fieldNums)) ):
    print('ERROR: Len field_DEC is '+str(len(field_DEC))+', len field_RA is '+str(len(field_RA))+', len field_fieldNum is '+str(len(field_fieldNums))+' and at least one of them doesn\'t match.')
    import sys
    sys.crash() #not a real call
#END IF
for i in range(0,field_numFields):
    dir_curr = os.path.join(dir_base,str(field_fieldNums[i]).zfill(2)) #create current path
    if os.path.isdir(dir_curr) == False:
        os.makedirs(dir_curr) #create dir if not a dir
    #END IF

    for j in range(0,field_numPerField):
        return_dict = astrobase.services.trilegal.query_radecl(field_RA[i], field_DEC[i], filtersystem='gaiaDR2', field_deg2=0.0841,
                    usebinaries=True, extinction_sigma=0, magnitude_limit=32.0, maglim_filtercol=1, trilegal_version=1.6,
                    extraparams=None, forcefetch=True, cachedir=dir_curr, verbose=True,
                    timeout=60.0, refresh=60.0, maxtimeout=1500.0) #read trilegal file, it is zipped

        with gzip.open(return_dict['tablefile'], 'rb') as f_in:
            with open(return_dict['tablefile'][:return_dict['tablefile'].find('.')]+str(j)+'.dat', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out) #un gzips and saves file, inspired by https://stackoverflow.com/a/44712152
            #END WITH
        #END WITH
    #END FOR j
    os.remove(return_dict['tablefile']) #cleanup (.gz file name is same for all FOR j runs, so only 1 .gz to cleanup at the end)
#END FOR i

print('done')