PRO nh_leisa_vega

  rad1 = 2.
  rad2 = rad1 * 2.
  
  data = MRDFITS('/data/NH/LEISA/lsb_0086619718_0x53c_sci.fit')

  lambda = MRDFITS('/data/NH/LEISA/lsb_0012474327_0x53d_sci_1.fit',1)
  lambda = REFORM(lambda[*,*,0])

  alpha_lyr = MRDFITS('~/Desktop/alpha_lyr_mod_002.fits',1)
  
  READCOL,'~/Desktop/leisa_vega_position.txt',frame,xpos,ypos,DELIMITER=','

  xpos = xpos - 1
  ypos = ypos - 1

  datsize = SIZE(data)
  xsize = datsize[1]
  ysize = datsize[2]
  nframes = datsize[3]

  ;data = data - REFORM(data[*,*,0])

  xframe = FLTARR(xsize,ysize)
  yframe = FLTARR(xsize,ysize)

  FOR xi=0,xsize-1 DO BEGIN
     xframe[xi,*] = xi
  ENDFOR
  
  FOR yi=0,ysize-1 DO BEGIN
     yframe[*,yi] = yi
  ENDFOR

  thisflux = !VALUES.F_NAN * FLTARR(nframes)
  thislam = !VALUES.F_NAN * FLTARR(nframes)
  psfwidth = !VALUES.F_NAN * FLTARR(nframes)
  
  FOR iframe=0,nframes-1 DO BEGIN

     IF FINITE(xpos[iframe]) THEN BEGIN
        rad = SQRT((xframe - xpos[iframe])^2 + (yframe - ypos[iframe])^2)
        whpl = WHERE(rad LE rad1,count1)
        ;thisdat = REFORM(data[xpos[iframe]-2:xpos[iframe]+2,$
        ;                      ypos[iframe]-2:ypos[iframe]+2,iframe])
        thisdat = REFORM(data[*,*,iframe] - data[*,*,0])
        whnt = WHERE(rad GT rad1 AND rad LE rad2,count2)
        thisflux[iframe] = (TOTAL(thisdat[whpl]) - $
                            TOTAL(thisdat[whnt[0:count1-1]]))
        ;ind = ARRAY_INDICES(thisdat,whpl)
        ;stop
        ;gf = GAUSS2DFIT(thisdat,A)
        ;psfwidth[iframe] = MEAN(A[2:3])
        thislam[iframe] = MEAN(lambda[whpl])
     ENDIF
     
  ENDFOR


  
  whpl = WHERE(FINITE(thisflux))
  thisflux = thisflux[whpl]
  thislam = thislam[whpl]
  ;whpl = WHERE(FINITE(psfwidth))
  ;print,median(psfwidth[whpl])
  
  thisomega = 1. * (12.55345/3600.)^2 * (!PI/180.)^2 ; Omega_pix
  thisomega1 = 1.13 * (18.1/3600.)^2 * (!PI/180.)^2  ; Omega_psf
  
  thisflux = thisflux * thisomega
  
  plot,thislam,thisflux
  oplot,alpha_lyr.wavelength/1e4,alpha_lyr.flux

  outspec = INTERPOL(alpha_lyr.flux,alpha_lyr.wavelength/1e4,thislam)

  oplot,thislam,outspec,col=34534545645

  plot,thislam,thisflux - outspec

  thisflux = thisflux * 1e-7 * 1e4 * 1e9 * (thislam * 1e4)

  outspec = outspec * 1e-7 * 1e4 * 1e9 * (thislam * 1e4)
  
  print,mean(thisflux - outspec)

  whpl = WHERE(thislam GT 1.8 AND thislam LT 2.4,count2)
  print,stddev(thisflux[whpl] - outspec[whpl]) / thisomega1 * $
        sqrt(0.586/4.) / FLOAT(count1)^(1.5);/SQRT(FLOAT(count2)); * 0.232

  ;; one sqrt factor of count1 comes from N_pix summed to get back per pixel rms
  ;; one factor of count1 comes from N_sum
  
  stop
  




END

