function [alpha,delta]=pix2radec(astrometry,xpix,ypix)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function pix2radec.m
%%  Sept. 14, 2010
%%  Mike Zemcov
%%
%%  This function takes a CIBER astrometry structure and a set pixel
%%  pairs and outputs the corresponding ra and dec of those pixels.  It
%%  follows the instructions here:
%%  http://ssc.spitzer.caltech.edu/files/spitzer/shupeADASS.pdf and here:
%%  
%%
%%  Inputs: astrometry - a CIBER-style astrometry structure
%%          xpix       - a list of x pixels
%%          ypix       - a list of y pixels
%%  Outputs: alpha     - ra corresponding to the input pixel vectors
%%           delta     - dec corresponding to the input pixel vectors
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  radeg = 180.0d0 ./ pi;
  
  u = xpix - astrometry.id__rpix1;
  v = ypix - astrometry.id__rpix2;

  up = u;  
  vp = v;
  
  if isfield(astrometry,'a_order') == 1 
  for ia=0:astrometry.a_order
    for ja=0:astrometry.a_order
      thisord = sprintf('a_%d_%d',ia,ja);
      if isfield(astrometry,thisord)
        up = up + getfield(astrometry,thisord) .* u.^(ia) .* v.^(ja);
      end
    end
  end
  
  for ia=0:astrometry.b_order
    for ja=0:astrometry.b_order
      thisord = sprintf('b_%d_%d',ia,ja);
      if isfield(astrometry,thisord)
        vp = vp + getfield(astrometry,thisord) .* u.^(ia) .* v.^(ja);
      end
    end
  end
  end
    
  %dbstop
  
  % Can't use matrix notation, in
  % case X and Y are vectors
  xsi = (astrometry.id_cd1_1 .* up + astrometry.id_cd1_2 .* vp);%.*(1-2.1e-3);
  eta = (astrometry.id_cd2_1 .* up + astrometry.id_cd2_2 .* vp);%.*(1-2.1e-3);

  wcsparams.lonpole = astrometry.id_lonpole;
  wcsparams.latpole = astrometry.id_latpole;
  wcsparams.crval = [astrometry.id_crval1,astrometry.id_crval2];
  
  [alpha,delta] = wcsxy2sph(xsi,eta,wcsparams);
  
  if numel(alpha) == astrometry.id_naxis1.*astrometry.id_naxis2
      alpha = reshape(alpha,astrometry.id_naxis2,astrometry.id_naxis1);
      delta = reshape(delta,astrometry.id_naxis2,astrometry.id_naxis1);
  end
    
  return

end
