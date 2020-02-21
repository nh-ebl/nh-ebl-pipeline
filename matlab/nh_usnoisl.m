function usnoisl = nh_usnoisl(data, paths)

  resize = 10;
   
  datafiles = dir(sprintf('%s*.mat',paths.datadir));
  
  load('lookup/nh_lorri_psf.mat');
    
  thispsf = psf.modelpsf;
  
  starimage = zeros(256.*resize);
  wingimage = zeros(256.*resize);
  
  load(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum));
    
  ncat = numel(B1mag);
    
  % loop over each catalog entry;
  for row = 1:ncat
  
    mags = [B1mag(row),B2mag(row),R1mag(row),R2mag(row),I2mag(row)];
    sg = [B1sg(row),B2sg(row),R1sg(row),R2sg(row),I2sg(row)];
    lambda_mag = [425,462.5,645,650,810];
    whpl = (mags < 21) & (mags > 1);
    if sum(whpl) > 1
      thismag = nh_synthetic_photometry(lambda_mag(whpl),mags(whpl),'LORRI');  
    else
      thismag = NaN;
    end
    
    %find x/y coordinate of the object
    [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
      
    if ypix >= 1 & ypix <= 256 & xpix >= 1 & xpix <= 256 & ~isnan(thismag)
      
      thisflux = data.cal.vzero .* 10^(-thismag./2.5);
      
      starimage(round(ypix.*resize),round(xpix.*resize)) = ...
	  starimage(round(ypix.*resize),round(xpix.*resize)) + thisflux;      
    
      if thismag < 17.75
	wingimage(round(ypix.*resize),round(xpix.*resize)) = ...
	    wingimage(round(ypix.*resize),round(xpix.*resize)) + thisflux;      
      end
    
    end
	
  end

  starimage_conv = conv2(starimage,thispsf,'same');
  wingimage_conv = conv2(wingimage,thispsf,'same');
    
  starimage_dec = resize.*imresize(starimage_conv,1./resize,'nearest');
  wingimage_dec = resize.*imresize(wingimage_conv,1./resize,'nearest');
    
  starimage_cal = starimage_dec .* data.cal.nu .* 1e-26 .* 1e9 ./ ...
      data.cal.omega_pix;
  wingimage_cal = wingimage_dec .* data.cal.nu .* 1e-26 .* 1e9 ./ ...
      data.cal.omega_pix;
  
  usnoisl.isltot = mean(starimage_cal(~data.mask.onemask));
  usnoisl.totimage = starimage_cal;
  usnoisl.islwing = mean(wingimage_cal(~data.mask.onemask));
  usnoisl.wingimage = wingimage_cal;
  usnoisl.islfaint = usnoisl.isltot - usnoisl.islwing;
  
%   figure(4); clf
%   imagesc(wingimage_cal.*~data.mask.onemask)
%   title(sprintf('%8.3f',usnoisl.islwing));
%   colorbar
%   caxis([0,10])
%   drawnow
    
end