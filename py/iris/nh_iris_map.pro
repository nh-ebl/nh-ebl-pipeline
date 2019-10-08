PRO nh_iris_map

  field = 5

  ;; path to the input files
  mappath = '../../../../IRISNOHOLES_B4H0'
  catpath = 'info_issa_map4.txt'

  ; set up the coordinates
  CASE field OF
     1: BEGIN
        ra  = 177.4586
        dec = 32.6273
     END
     2: BEGIN
        ra = 196.0334
        dec = 23.9452
     END
     3: BEGIN
        ra  = 161.8978
        dec = -26.7821
     END
     4: BEGIN
        ra = 177.5121
        dec = 32.6308
     END
     5: BEGIN
        ra = 196.0075
        dec = 23.9506
     END
     6: BEGIN
        ra = 346.1128
        dec = -7.166
     END
     7: BEGIN
        ra = 1.8064
        dec = -1.2497
     END
     8: BEGIN
        ra = 270.6523
        dec = -14.6237
     END
     9: BEGIN
        ra = 237.5430
        dec = -33.8045
     END
     10: BEGIN
        ra = 268.4866
        dec = -34.7916
     END
     11: BEGIN
        ra = 269.6507
        dec = -15.4847
     END
     12: BEGIN
        ra = 231.3464
        dec = -17.9611
     END
     13: BEGIN
        ra = 269.7189
        dec = -15.4893
     END
     14: BEGIN
        ra = 237.3581
        dec = -33.7934
     END
  ENDCASE

  ;; default pixel size of LORRI
  pixsize = 4.1 / 3600.

  ;; want a patch as large as the maximum possible extent of the image
  naxis = FIX(SQRT(2.) * [256,256])

  ;; this is a dummy holder for the x,y coords we are interested in
  xvec = REPVEC(FINDGEN(naxis[0]),NLEN=naxis[1])
  yvec = TRANSPOSE(REPVEC(FINDGEN(naxis[1]),NLEN=naxis[0]))

  ;; make the astrometry structure
  MAKE_ASTR,astr,CD=[[-pixsize,0.],[0.,pixsize]],CRPIX=naxis / 2 + 1,$
            CRVAL=[ra,dec]

  ;; now make the vectors of the longitude and latitude
  XY2AD,xvec,yvec,astr,lonvec,latvec

  ;; make a dummy fits header to put info into
  MKHDR,header,latvec

  ;; put the header info into the dummy
  PUTAST,header,astr

  ;; the next step takes some time so warn me that I've
  ;; actually done something good
  MESSAGE,'Generated astrometry info, getting IRAS values.',/INFORMATIONAL

  ;; get the IRAS values at the astrometry positions
  map = MOSAIQUE_IRIS(header,band=4,catname=catpath,dir=mappath)

  ;; and write out the output
  ; commented out this part while I do test things
  ; WRITEFITS,'/data/zemcov/NH/IRIS/iris_'+STRING(field,FORMAT="(I02)")+'_fx.fits',$
  ;           map,header
  ; WRITEFITS,'/data/zemcov/NH/IRIS/iris_'+STRING(field,FORMAT="(I02)")+'_ra.fits',$
  ;           lonvec,header
  ; WRITEFITS,'/data/zemcov/NH/IRIS/iris_'+STRING(field,FORMAT="(I02)")+'_dc.fits',$
  ;           latvec,header

END
