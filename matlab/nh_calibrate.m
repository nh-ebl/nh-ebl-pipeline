%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_calibrate.m
%%  Jun 2016, MZ
%%  This function applies known calibration factors to the reducing NH data.
%%  Input = data, an NH pipeline data structrue.
%%  Output = data, an NH pipeline data structure.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = nh_calibrate(data)

  % some constants
  data.const.Jy = 1e-26; % W m^-2 sr^-1
  data.const.nW = 1e9; % W per nW

  % omega of beam compared to one pixel
  data.cal.pixperbeam = 2.640;
  
  % multiplicative zero point of the conversion
  data.cal.szero = 10.^(data.cal.magoff ./ -2.5);
  
  % zero point of the vega system in R_L band
  data.cal.vzero = 3050; % Jy at R_L=0.
  
  % compute the number of Jy/bit = 10^-26 W m^-2 Hz^-1 / bit
  data.cal.Jyperbit = data.cal.vzero .* data.cal.szero; % Jy/bit

  % some astrometric quantities we'll need that are fixed values from elsewhere
  data.cal.pixsize = 1.05.*4.*4.96e-6; % radians
  data.cal.pixsize_arcsec = data.cal.pixsize .* 180 .* 3600 ./ pi; % arcsec
  data.cal.omega = data.cal.pixperbeam .* data.cal.pixsize.^2;
  data.cal.omega_pix = (data.cal.pixsize_arcsec./3600.).^2 * (pi./180.).^2;
  
  % compute the frequency in Hz
  data.cal.nu = 3e8 ./ (data.header.lambda .* 1e-9);
  
  % so this becomes the conversion from DN/s to nW m^-2 sr^-1
  data.cal.sbconv = data.cal.Jyperbit * data.const.Jy .* data.const.nW .* ...
      data.cal.nu ./ data.cal.omega;
  
  load('lookup/nh_refcorr.mat');
  
  whpl = data.header.date_jd == correction.date;

  % calculate some of the stats to take account of the calibration factor
  data.stats.calmean = data.cal.sbconv .* data.stats.maskmean;
  data.stats.calerr = data.cal.sbconv .* data.stats.maskerr;
  
  data.stats.corrmean = data.cal.sbconv .* (data.stats.maskmean + ...
      correction.corr(whpl)./data.header.exptime);
  data.stats.correrr = data.cal.sbconv .* (data.stats.maskerr);
    
  % pull out the raw image data
  rawimage = data.data;
    
  % and make calibrated versions of that as well.
  data.image.rawimage = rawimage ./ data.header.exptime;
  data.image.calimage = rawimage .* data.cal.sbconv ./ data.header.exptime;

% end


