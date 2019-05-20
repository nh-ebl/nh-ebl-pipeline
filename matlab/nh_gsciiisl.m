function gsciiisl = nh_gsciiisl(data)

  resize = 10;

  paths = get_paths();
   
  datafiles = dir(sprintf('%s*.mat',paths.datadir));
  
  load('lookup/nh_lorri_psf.mat');
  
  thispsf = imresize(psf.psf,resize,'nearest');
  
  thispsf = psf.modelpsf;
  
  starimage = zeros(256.*resize);
  wingimage = zeros(256.*resize);

  filename = sprintf('%sgscii_%02dred.txt.csv',...
      paths.gscdir,data.header.fieldnum);

  [ra,dec,UMag,BMag,VMag,RMag,IMag] = ...
      textread(filename,'%f,%f,%f,%f,%f,%f,%f');

  whpl = VMag < 30;
  
  ra = ra(whpl);
  dec = dec(whpl);
  VMag = VMag(whpl);
      
  ncat = numel(VMag);
    
  % loop over each catalog entry;
  for row = 1:ncat
  
    thismag = VMag(row);
    
    %find x/y coordinate of the object
    [ypix, xpix] = radec2pix(ra(row),dec(row), data.astrom);
      
    if ypix >= 1 & ypix <= 256 & xpix >= 1 & xpix <= 256 & ~isnan(thismag)
      
      thisflux = 3636 .* 10^(-thismag./2.5);
      
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
    
  gsciiisl.isltot = mean(starimage_cal(~data.mask.onemask));
  gsciiisl.totimage = starimage_cal;
  gsciiisl.islwing = mean(wingimage_cal(~data.mask.onemask));
  gsciiisl.wingimage = wingimage_cal;
  gsciiisl.islfaint = gsciiisl.isltot - gsciiisl.islwing;
    
  %figure(1); clf
  %imagesc(wingimage_cal.*~data.mask.onemask)
  %title(sprintf('%8.3f',gsciiisl.islwing));
  %colorbar
  %caxis([0,10])
  %drawnow
  
  
% end