function [Z] = make_maskslow( filename, band, catalog )
% Create star mask of input image using input catalog
%       Loops through each pixel of image, determines RA/DEC of pixel.
%       Checks if there is a object at the location by searching though
%       catalog. If there is determines radius of object and adds object to
%       mask of the corresponding counts.
%
% Old method - use make_mask.m.
%
% PARAMETERS:
%       filename = .fits file
%       band = passband of mask
%       catalog = .mat file containing arrays for catalog data
%                 [RA, DEC, B1mag, B2mag, R1mag, R2mag]
%
% OUTPUT:
%        Z = an array mask for filename. This will can be removed since file
%	     is written.
%
% author: Poppy Immel
% date: May 3, 2016
% email: pgi8114@rit.edu

%Initialize mask array. Size of mask: 256x256 pix. Corresponds to
%image size
xdim = 256;
ydim = 256;
X = 1:xdim;
Y = 1:ydim;
Z = fitsread(filename);%zeros(ydim,xdim);

%Load cooresposing catalog matrix file, inlcudes arrays [RA, DEC, B1mag,
%B2mag, R1mag, R2mag]
load(catalog)

%Set catalog ra/dec arrays and determine size of catalog
CAT_RA = RA;
CAT_DEC = DEC;
[~,ncat] = size(RA);

%create astromety type object to use for pix2radec funtion.
info = fitsinfo(filename);
astrm = fits_to_astrometry(info);

%calculate pixel size in degrees
%may need to ajust since the RA/Dec axis are not the same as the x/y axis
[RA11,DEC11] = pix2radec(astrm, 1, 1);
[RA12,~] = pix2radec(astrm, 1, 2);
[~,DEC21] = pix2radec(astrm, 2, 1);
RArange = abs(RA11-RA12)/2;
DECrange = abs(DEC11-DEC21)/2;

%Set mag as magnatude array corresponing to the band specified
mag = B1mag;
if strcmp(band, 'B2')
    mag = B2mag;
elseif strcmp(band, 'R1')
    mag = R1mag;
elseif strcmp(band, 'R2')
    mag = R2mag;
end

%Set values for radius calculation: r(mag) = -alpha(mag - max_mag)+min_r
min_r = 1;
max_mag = 16; %needs to be changed to limiting mag of instrument.
alpha = 1; %needs to be calculated

%loop over each x/y pixel
for xpix = X
    for ypix = Y
        %convert pix to ra/dec
        [currRA,currDEC] = pix2radec(astrm, xpix, ypix);
        %search to ra/dec in catalog
        for colcat = 1:ncat
            %check if mag for current band is 0, if it is, ignore
            if mag(colcat) == 0
                continue;
            end
            %find objects ra/dec range for current pixel.
            if currRA >= CAT_RA(colcat)-RArange ...
                    && currDEC <= CAT_DEC(colcat)+DECrange ...
                    && currDEC >= CAT_DEC(colcat)-DECrange ...
                    && currRA <= CAT_RA(colcat)+RArange
                
                %set magnitude
                curr_mag = mag(colcat);
                
                %determine radius and ints(counts) of object
                radius = round(-alpha*(curr_mag - max_mag) + min_r);
                %ints = mag_to_int(curr_mag, band);
                
                %create submask of object (create array of 0's and 1's,
                %where the 1's represnt to location of the objects in the
                %submask). Basically a sqaure of 0's with a circle of 1's
                [rr, cc] = meshgrid(1:2*radius+1);
                C = sqrt((rr-radius-1).^2+(cc-radius-1).^2)<=radius;
                
                %combined the submask (C) and mask (Z) where xpix, ypix is
                %the center of the object
                for i = 1:(2*radius+1)
                    xcurr = xpix-radius-1+i;
                    if xcurr < 1 || xcurr > xdim;
                        continue;
                    end
                    for j = 1:(2*radius+1)
                        ycurr = ypix-radius-1+j;
                        if ycurr < 1 || ycurr > ydim;
                            continue;
                        end
                        %check that this is the
                        %brightest object in this
                        %pixel.
                        %if C(i,j) == 0;
                        %	continue;
                        %end
                        %if ints > Z(ycurr,xcurr)
                        %	Z(ycurr,xcurr) = ints*C(i,j);
                        if C(i,j) ~= 0
                            Z(ycurr,xcurr) = NaN;
                        end
                    end
                end
            end
        end
    end
end

%%write masked image
% newfits = sprintf('%s_masksubtracted.fits', band);
% imagedata = single(Z);
% fitswrite(imagedata, newfits);
end


