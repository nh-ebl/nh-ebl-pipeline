function nh_calc_maskfluxfrac()

  pdim = 1001;
  resize = 10;
  mags = [5:20];
  nmags = numel(mags);
  max_mag = 17.75;
  
  sigma = resize .* 1.5 ./ (2.*sqrt(2.*log(2)))
  
  y = gauss2([1,(pdim-1)/2-4.5,(pdim-1)/2-4.5,sigma],[1:pdim],[1:pdim]');
  y = y ./ sum(y(:));
  
  [xp,yp] = meshgrid([-50:50],[-50:50]);
  rad = sqrt(xp.^2 + yp.^2);

  for imag=1:nmags
    %thisimage = zeros(pdim);
    thisimage = 3050.*10^(-mags(imag)./2.5).*y;
  
    thisimage = imresize(thisimage,1./resize);
    
    figure(1); clf
    imagesc(thisimage)
    colorbar
    %xlim([40,60])
    %ylim([40,60])
    drawnow
    pause(0.5)
  
    %radius = 6000.*exp(-mags(imag))+2.5;
    radius = 2.5*(max_mag./mags(imag)).^(2);%
    whpl = rad < radius;
    thisimagep = thisimage;
    thisimagep(whpl) = 0;
    
    figure(1); clf
    imagesc(thisimagep)
    colorbar
    %xlim([40,60])
    %ylim([40,60])
    drawnow
    pause(0.5)
    
    mags(imag)
    100.*abs(sum(thisimagep(:))./sum(thisimage(:)))
    abs(mean(thisimage(~whpl))).*(3e8./655e-9) .* 1e-26 .* 1e9 ./ ((4.3./3600).^2.*(pi./180).^2)
    
    
    
  
  end
    
  dbstop
    







% end