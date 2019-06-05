%This function performs background-subtracted aperture photometry
%Inputs: image containing source, source pixel coordinates, desired radius
%of aperture, and desired multiplicative factor for sky annulus (ninner*rad
%= inner annulus rad, nouter*rad = outer annulus rad)

%Symons June 2019

function skysub_apsum = ap_photom(image, sourcex, sourcey, rad, ninner, nouter)

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
skysub_apsum = apsum - bkgsum;


%inputs like data.data.*~data.mask.onemask, data.ghost.ghostx,
%data.ghost.ghosty, data.ghost.ghostrad
%ninner = 2 nouter = 3
end