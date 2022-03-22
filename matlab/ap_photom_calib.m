%This function performs background-subtracted aperture photometry
%Inputs: image containing source, source pixel coordinates, desired radius
%of aperture, and desired multiplicative factor for sky annulus (ninner*rad
%= inner annulus rad, nouter*rad = outer annulus rad)

%inputs like data.data.*~data.mask.onemask, data.ghost.ghostx,
%data.ghost.ghosty, data.ghost.ghostrad
%ninner = 2 nouter = 3

%Symons June 2019

function skysub_apsum = ap_photom_calib(star_image, bkg_image, sourcex, sourcey, rad, ninner, nouter)

%Aperture photometry method to determine background-subtracted counts for
%ghost inside aperture
%meshgrid based on image size
[xgrid, ygrid] = meshgrid(1:size(star_image,2), 1:size(star_image,1));
%create mask where meshgrid values are less than radius distance from
%source
mask = ((xgrid-sourcex).^2 + (ygrid-sourcey).^2) <= rad.^2;
%list of all values in masked area
values = star_image(mask);
%aperture sum of all values inside aperture
apsum = sum(values);

%now calculating sky background
[xgrid, ygrid] = meshgrid(1:size(star_image,2), 1:size(star_image,1));
%create mask where values are within outer sky annulus radius
skymaskinner = ((xgrid-sourcex).^2 + (ygrid-sourcey).^2) <= (rad*nouter).^2;
%create mask where values are outside inner sky annulus radius
skymaskouter = ((xgrid-sourcex).^2 + (ygrid-sourcey).^2) >= (rad*ninner).^2;
%combine intersection of those masks
skymask = skymaskinner & skymaskouter;
%list of values in masked area
skyvalues = bkg_image(skymask);
%sum of all values in sky annulus
skysum = sum(skyvalues);

%calculate area of sky annulus
skymaskarea = pi*((rad*nouter).^2 - (rad*ninner).^2);
%calculate area of aperture
aparea = pi*rad.^2;

%compute background summed value scaled for difference in areas
bkgsum = skysum/skymaskarea*aparea;
%compute the background subtracted aperture sum
skysub_apsum = apsum - bkgsum;
% flux = skysub_apsum/data.astrom.exptime;
% mag = -2.5*log10(flux)+20;
% info = '';
% info = [info,'Background-subtracted counts: ',num2str(skysub_apsum),' Flux: ',num2str(flux),' Mag: ',num2str(mag)];


%Plot entire image with circle around ghost
h = figure(1);
clf;
% set(h,'visible','off');
th = 0:pi/50:2*pi;
xunit = rad*cos(th)+sourcex;
yunit = rad*sin(th)+sourcey;
imagesc(star_image);
hold on;
plot(xunit,yunit,'r');
pbaspect([1 1 1]);
a = colorbar;
a.Label.String = 'Intensity [DN]';
caxis([0,5]);
% title(sprintf('Ghost Mean: %.2f Bkg: %.2f Sub Mean: %.2f Num Pix: %d Ghost Tot: %.2f',mean(image(mask & image~=0)),bkg,pixval,length(values),skysub_apsum));
set(gca,'YDir','normal');
xlabel('LORRI X Pixels');
ylabel('LORRI Y Pixels');
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.magdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');
end