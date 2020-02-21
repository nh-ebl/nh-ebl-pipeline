function data = nh_calcdgl(data, paths)
  
  dglparams = nh_get_dgl_params();
  
  irismap = fitsread(sprintf('%siris_%02d_fx.fits',...
      paths.irisdir,data.header.fieldnum));
  irisra = fitsread(sprintf('%siris_%02d_ra.fits',...
      paths.irisdir,data.header.fieldnum));
  irisdc = fitsread(sprintf('%siris_%02d_dc.fits',...
      paths.irisdir,data.header.fieldnum));
  
  %create interpolation of iris map and lorri field
  F = TriScatteredInterp(irisra(:),irisdc(:),irismap(:));
  irisim = F(data.astrometry.ra,data.astrometry.dec);
  
  %calculte mean and std of iris map
  ohm_mean = nanmean(irisim(:));
  ohm_std = nanstd(irisim(:));
  
  %[lorri_lambda,lorri_trans] = textread('lookup/nh_lorri_bandpass.txt','%f %f');
  
  %[dgl_lambda,dgl_conv] = textread('lookup/nh_dgl.txt','%f, %f');
  
  %dgl_lambda = 1e3 .* dgl_lambda;
  %dgl_conv = dgl_conv ./ 2;
  
  %lambda = [300:1000];
  
  %myconv = interp1(dgl_lambda,dgl_conv,lambda,'spline');
    
  %myconv = sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0) .* ...
  %    myconv) ./ ...
  %    sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0));
  
  %myfreq = 3e8 ./ 100e-6;
  
  %nu I_nu is (iris map - 0.8 MJy/sr) * 1e-20 * 3e8/100e-6 * 1e9
  nuinu = (irisim-dglparams.cib(1)) .* dglparams.norm;
  
  
  %b_lambda (I_nu_opt/I_nu_100um) is estimated by fitting mean ZDA04 model to many measurements
  %of b_lambda
  %c_lambda (lambda*I_lambda_opt/nu*I_nu_100um)s is 10^-6*(100um/0.655um)*b_lambda
  %cbar is 0.491 (cbar_655nm calculated by integrating cbar over lorri
  %bandpass)
  cbar = dglparams.cbar(1);
  
  %A = 1/0.567 is normalizing factor d0 (from normalizing d(b) at b=25)
  %g = 0.61 is bandpass-weighted mean of observation-constrained model for
  %high-lat diffuse dust component of DGL
  %dl is d(b), a function that accounts for change in c_lambda due to b
  dl = dglparams.A .* (1 - 1.1 .* dglparams.g(1).*sqrt(sin(abs(data.coords.galactic(2)).*pi./180)));
  
  %dglim is estimate of DGL using average of 100um intensity, c_lambda, and
  %d(b)
  dglim = nuinu .* cbar .* dl;
  
  %propagated error here
  dglerr_nuinu = (cbar .* dl).^2 .* dglparams.cib(2).^2;
  dglerr_cbar = (nanmean(nuinu(:)) .* dl).^2 .* dglparams.cbar(2).^2;
  dglerr_dl = (nanmean(nuinu(:)) .* cbar .* dglparams.A .* 1.1 .* ...
      sqrt(sin(abs(data.coords.galactic(2)).*pi./180))).^2 .* dglparams.g(2).^2;
  
  dgl_err = sqrt(dglerr_nuinu + dglerr_cbar + dglerr_dl);
  
  %whpl = irisim > 10;
  %load('lookup/nh_dglcorr.mat');
  %corrfac = interp1(dglcorr.ohm,dglcorr.corr,irisim(whpl));
  %dglim(whpl) = corrfac .* dglim(whpl);
  
  %mean and std of dgl estimate (at each pixel?)
  dgl_mean = nanmean(dglim(:));
  dgl_std = nanstd(dglim(:));
  
  data.dgl.dglmean = dgl_mean;
  data.dgl.dglstd = dgl_std;
  data.dgl.dglerr = dgl_err;
  data.dgl.ohmmean = ohm_mean;
  data.dgl.ohmstd = ohm_std;
  data.dgl.ohmim = irisim;
  data.dgl.dglim = dglim;
  data.dgl.conv = dglparams.cbar(1).*dl;
  data.dgl.convp = dglparams.norm .* dglparams.cbar(1).*dl;
    
  %disp(sprintf('Field number: %d; b: %4.2f; DGL: %7.3f; DGLerr: %7.3f.',...
  %    data.header.fieldnum,data.coords.galactic(2),dgl_mean,dgl_err))
  
% end