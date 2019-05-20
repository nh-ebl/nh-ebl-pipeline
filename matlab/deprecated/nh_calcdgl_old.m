function data = nh_calcdgl(data)

  paths = get_paths();
  
  irismap = fitsread(sprintf('%siris_%02d_fx.fits',...
      paths.irisdir,data.header.fieldnum));
  irisra = fitsread(sprintf('%siris_%02d_ra.fits',...
      paths.irisdir,data.header.fieldnum));
  irisdc = fitsread(sprintf('%siris_%02d_dc.fits',...
      paths.irisdir,data.header.fieldnum));
  
  F = TriScatteredInterp(irisra(:),irisdc(:),irismap(:));
  irisim = F(data.astrometry.ra,data.astrometry.dec);

  ohm_mean = nanmean(irisim(:));
  ohm_std = nanstd(irisim(:));
  
  [lorri_lambda,lorri_trans] = textread('lookup/nh_lorri_bandpass.txt','%f %f');
  
  [dgl_lambda,dgl_conv] = textread('lookup/nh_dgl.txt','%f, %f');
  
  dgl_lambda = 1e3 .* dgl_lambda;
  %dgl_conv = dgl_conv ./ 2;
  
  lambda = [300:1000];
  
  myconv = interp1(dgl_lambda,dgl_conv,lambda,'spline');
  
  whpl = abs(lambda-640) == min(abs(lambda-640));
  
  % factor 0.357 comes from Brandt & Draine
  myconv = myconv ./ myconv(whpl) .* 0.357;
    
  myconv = sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0) .* ...
      myconv) ./ ...
      sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0));
  
  myfreq = 3e8 ./ 100e-6;
  
  dglim = (irisim-0.8) .* 1e-20 .* myfreq .* myconv .* 1e9;
  
  %whpl = irisim > 10;
  %load('lookup/nh_dglcorr.mat');
  %corrfac = interp1(dglcorr.ohm,dglcorr.corr,irisim(whpl));
  %dglim(whpl) = corrfac .* dglim(whpl);
  
  dgl_mean = nanmean(dglim(:));
  dgl_std = nanstd(dglim(:));
  
  data.dgl.dglmean = dgl_mean;
  data.dgl.dglstd = dgl_std;
  data.dgl.ohmmean = ohm_mean;
  data.dgl.ohmstd = ohm_std;
  data.dgl.ohmim = irisim;
  data.dgl.dglim = dglim;
  data.dgl.conv = myconv;
  data.dgl.convp = 1e-20 .* myfreq .* myconv .* 1e9;
    
% end