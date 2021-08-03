%This function takes a data (.mat) file for one LORRI image and calculates
%the diffuse contribution from optical ghosting from all stars within range
%to cause a ghost. This list was previously determined in nh_findghoststar.

%The contribution from all stars is added together to get the total diffuse
%ghost contribution for that image. This is later subtracted from the image
%mean.

%Symons 2021

function [data] = nh_diffghost(data,paths)

% Calculate number of potential bright stars contributing to ghost
numstars = size(data.ghost.brightmag,1);

% Preallocate space for list of star distance to center and x and y coord of star
stardistcent = zeros(1,numstars);

% For all stars, retrieve x and y position
for j = 1:numstars
    
    x = data.ghost.brightxpix(j,1);
    y = data.ghost.brightypix(j,1);
    
    % Adjust coordinates for large (654x654) grid
    
    x = 199 + x;
    
    y = 199 + y;
    
    % Calculate distance from center pixel to star pixel and
    % distance from star pixel to ghost pixel
    stardistcent(1,j) = sqrt((x-(199+128))^2 + (y-(199+128))^2);
end

% Find which star (if more than one) is closest to center (assuming
% that is the cause of the ghost)
[M,I] = min(stardistcent);

% For all stars, retrieve x and y position (except star that is causing
% masked ghost)
boxxpix = data.ghost.boxxpix(data.ghost.boxxpix~=data.ghost.brightxpix(I,1)) ;
boxypix = data.ghost.boxypix(data.ghost.boxypix~=data.ghost.brightypix(I,1));
boxmag = data.ghost.boxmag(data.ghost.boxmag~=data.ghost.brightmag(I,1));

% Adjust coordinates for large (654x654) grid
xbox = 199 + boxxpix;

ybox = 199 + boxypix;

% Calculate distance from center pixel to star pixel and
% distance from star pixel to ghost pixel
stardistcentbox = sqrt((xbox-(199+128)).^2 + (ybox-(199+128)).^2);

% Cut down list to stars in range to cause ghost
boxxpixcut = boxxpix(stardistcentbox <= 268);
boxypixcut = boxypix(stardistcentbox <= 268);
boxmagcut = boxmag(stardistcentbox <= 268);

% Calculate number of pixels in masked ghost area and number of pixels in
% image
[xgrid, ygrid] = meshgrid(1:size(data.data,2), 1:size(data.data,1));
maskboxtest = ((xgrid-(100+128)).^2 + (ygrid-(100+128)).^2) <= 21.5.^2;
ghost_pix = length(data.data(maskboxtest));
image_pix = data.astrom.imagew.*data.astrom.imageh;

% Calculate total surface brightness of all diffuse ghosts from all stars
ghostsb = (10.^(-0.3248.*boxmagcut + 6.4991))./image_pix; % These values determined in ghost_analysis.m from fit
ghostsberrpos = ((10.^((-0.3248.*boxmagcut + 6.4991)+0.0733))./image_pix) - ((10.^((-0.3248.*boxmagcut + 6.4991)))./image_pix); % These values determined in ghost_analysis.m from fit
ghostsberrneg = ((10.^((-0.3248.*boxmagcut + 6.4991)))./image_pix) - ((10.^((-0.3248.*boxmagcut + 6.4991)-0.0733))./image_pix); % These values determined in ghost_analysis.m from fit

% Calculate surface brightness from stars that cause ghosts
Fcat = (data.cal.vzero) .* 10.^(-(boxmagcut+data.cal.gaia2lorrimag)./2.5);
surveyarea = data.astrom.imagew.*data.astrom.imageh.*...
    data.cal.pixsize_arcsec.^2 ./ 3600.^2;
starsb = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea .* (pi./180).^2);

% Save summed diffuse ghost surface brightness to data
data.ghost.diffusesub = sum(ghostsb);
data.ghost.diffusesuberrpos = sqrt(sum(ghostsberrpos.^2));
data.ghost.diffusesuberrneg = sqrt(sum(ghostsberrneg.^2));

% Calculate slope of ghost sb vs. star sb
fitter = linear_fit(log10(starsb),log10(ghostsb));
starfit=linspace(min(log10(starsb)),max(log10(starsb)));
fitterlin = linear_fit(10.^(starfit)',10.^(fitter.m*starfit' + fitter.b));
ghostfit=(fitter.m*starfit + fitter.b);

% Save slope to data to compare to ISL
data.ghost.diffuseslope = fitterlin.m;

% Plot ghost sb vs. star sb with linear fit
% h = figure(1);
% clf;
% set(h,'visible','off');
% scatter(starsb,ghostsb,'filled');
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% hold on;
% fit = plot(10.^(starfit),10.^(ghostfit));
% ylabel('Ghost Intensity [nW m^{-2} sr^{-1}]');
% xlabel('Star Intensity [nW m^{-2} sr^{-1}]');
% title(sprintf('Linear fit: y = %.5fx + %.5f',fitterlin.m,fitterlin.b));
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.ghostdiffdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

% t = tiledlayout(1,1);
% ax1 = axes(t);
% scatter(ax1,log10(starsb),log10(ghostsb),'r')
% ax1.XColor = 'r';
% ax1.YColor = 'r';
%
% ax2 = axes(t);
% scatter(ax2,boxmagcut,log10(ghostsb),'k')
% set(gca, 'XDir','reverse')
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% ax1.Box = 'off';
% ax2.Box = 'off';

end