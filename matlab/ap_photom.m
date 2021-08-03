%This function performs background-subtracted aperture photometry
%Inputs: image containing source, source pixel coordinates, desired radius
%of aperture, and desired multiplicative factor for sky annulus (ninner*rad
%= inner annulus rad, nouter*rad = outer annulus rad)

%inputs like data.data.*~data.mask.onemask, data.ghost.ghostx,
%data.ghost.ghosty, data.ghost.ghostrad
%ninner = 2 nouter = 3

%Symons June 2019

function skysub_apsum = ap_photom(image, sourcex, sourcey, rad, ninner, nouter, data, paths)

%Aperture photometry method to determine background-subtracted counts for
%ghost inside aperture
%meshgrid based on image size
[xgrid, ygrid] = meshgrid(1:size(image,2), 1:size(image,1));
%create mask where meshgrid values are less than radius distance from
%source
mask = ((xgrid-sourcex).^2 + (ygrid-sourcey).^2) <= rad.^2;
%list of all values in masked area
values = image(mask);
%aperture sum of all values inside aperture
apsum = sum(values);

%now calculating sky background
[xgrid, ygrid] = meshgrid(1:size(image,2), 1:size(image,1));
%create mask where values are within outer sky annulus radius
skymaskinner = ((xgrid-sourcex).^2 + (ygrid-sourcey).^2) <= (rad*nouter).^2;
%create mask where values are outside inner sky annulus radius
skymaskouter = ((xgrid-sourcex).^2 + (ygrid-sourcey).^2) >= (rad*ninner).^2;
%combine intersection of those masks
skymask = skymaskinner & skymaskouter;
%list of values in masked area
skyvalues = image(skymask);
%sum of all values in sky annulus
skysum = sum(skyvalues);

%calculate area of sky annulus
skymaskarea = pi*((rad*nouter).^2 - (rad*ninner).^2);
%calculate area of aperture
aparea = pi*rad.^2;

%compute background summed value scaled for difference in areas
bkgsum = skysum/skymaskarea*aparea;
%compute the background subtracted aperture sum
% skysub_apsum = apsum - bkgsum;
% flux = skysub_apsum/data.astrom.exptime;
% mag = -2.5*log10(flux)+20;
% info = '';
% info = [info,'Background-subtracted counts: ',num2str(skysub_apsum),' Flux: ',num2str(flux),' Mag: ',num2str(mag)];

% Make ghost histogram as alternative method of estimating number of counts
% inside ghost
%Preallocate for edge values of peak histogram value
minmaxedge = zeros(2,1);
%Calculate background (mean of unmasked pixels, including masked ghost)
bkg = mean(data.image.calimage(~data.mask.onemask));
%Save only ghost area of image
idx = image~=0;
ghost = image(mask&idx)-bkg;
%Find indices where pixel values are not 0 (masked)
idx = ghost~=0;
%Save the bin edges and bin counts for a histogram of the ghost
[N,edges] = histcounts(ghost(idx));
%Find the bin with maximum counts
[M,I] = max(N);
%Calculate edge values for that bin
minmaxedge(1,1) = edges(I);
minmaxedge(2,1) = edges(I+1);
%Pixel value is average of those bin edges
% pixval = median(minmaxedge);
% Mean ghost pixel value - hist seems Poisson-like
pixval = mean(ghost(idx));
mean(minmaxedge)-median(minmaxedge);
if pixval == 0.5
    fprintf('Ghost is not completely masked')
end
%Total counts for ghost are pixel value x number of pixels in ghost
skysub_apsum = pixval*length(values);
flux = skysub_apsum/data.astrom.exptime;
mag = -2.5*log10(flux)+20;
info = '';
info = [info,'Background-subtracted counts: ',num2str(skysub_apsum),' Flux: ',num2str(flux),' Mag: ',num2str(mag)];

%Plot actual histogram of ghost
% h = figure(1);
% clf;
% set(h,'visible','off');
% g = histogram(ghost(idx));
% hold on;
% xl1 = xline(median(ghost(idx)),'r');
% xl1.LineWidth = 2;
% xl2 = xline(mean(ghost(idx)),'b');
% xl2.LineWidth = 2;
% legend([xl1,xl2],{'Median','Mean'},'Location','northeast');
% title(sprintf('%s',data.header.rawfile));
% xlabel('Intensity [nW m^{-2} sr^{-1}]');
% ylabel('N');
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.ghosthistdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

%Plot 'bullseye' of ghost aperture and sky annulus
% h = figure(1);
% clf;
% set(h,'visible','off');
% totalmask = skymask + mask;
% imagesc(totalmask.*image);
% pbaspect([1 1 1]);
% colorbar;
% caxis([-20,20]);
% % t2 = annotation('textbox',[0.13,0.075,0,0],'string',info,'FitBoxToText','on');
% set (gcf, 'WindowButtonMotionFcn', @mouseMove);
% colorbar; 
% caxis([-20,20]);
% grid on;
% 
% t2 = annotation('textbox',[0.13,0.075,0,0],'string',info,'FitBoxToText','on'); 
% t2.LineStyle = 'none';
% title(sprintf('%s',data.header.rawfile));
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.magdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

%Plot entire image with circle around ghost
% h = figure(2);
% clf;
% set(h,'visible','off');
% th = 0:pi/50:2*pi;
% xunit = rad*cos(th)+sourcex;
% yunit = rad*sin(th)+sourcey;
% imagesc(image);
% hold on;
% plot(xunit,yunit);
% pbaspect([1 1 1]);
% a = colorbar;
% a.Label.String = 'Intensity [nW m^{-2} sr^{-1}]';
% caxis([-250,250]);
% % t2 = annotation('textbox',[0.13,0.075,0,0],'string',info,'FitBoxToText','on');
% % t2.LineStyle = 'none';
% % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
% % colorbar; 
% % grid on;
% % title(sprintf('%s',data.header.rawfile));
% title(sprintf('Ghost Mean: %.2f Bkg: %.2f Sub Mean: %.2f Num Pix: %d Ghost Tot: %.2f',mean(image(mask & image~=0)),bkg,pixval,length(values),skysub_apsum));
% set(gca,'YDir','normal');
% xlabel('LORRI X Pixels');
% ylabel('LORRI Y Pixels');
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.magdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

end