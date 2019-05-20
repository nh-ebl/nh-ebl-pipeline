function [Z] = make_mask( filename, band, catalog )
% Create star mask of input image using input catalog
%       Loops through each pixel of image, determines RA/DEC of pixel.
%       Checks if there is a object at the location by searching though
%       catalog. If there is determines radius of object and adds object to
%       mask of the corresponding counts.
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

%Load fits file and store as a matrix
Z = fitsread(filename);
[xdim, ydim] = size(Z);

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

%Set mag as magnatude array corresponing to the band specified
MAG = B1mag;
if strcmp(band, 'B2')
    MAG = B2mag;
elseif strcmp(band, 'R1')
    MAG = R1mag;
elseif strcmp(band, 'R2')
    MAG = R2mag;
end

%Set values for radius calculation: r(mag) = -alpha(mag - max_mag)+min_r
min_r = 1;
max_mag = 18; %needs to be changed to limiting mag of instrument.
alpha = .75; %needs to be calculated

%look over each catalog entry;
for row = 1:ncat
    %find x/y coordinate of ech object in catalog
    
    %TEMP VALUES FOR TESTING
    % xpix = 10; ypix = 10;
    [ypix, xpix] = radec2pix(CAT_RA(row),CAT_DEC(row), astrm);
                        %FUNCTION REQUIRES wcssph2xy.m
    %check if the object is in the image
    if xpix > 0 && xpix <= xdim && ypix > 0 && ypix <= ydim
        xpix = round(xpix);
        ypix = round(ypix);
        
        curr_mag = MAG(row);
        if curr_mag > max_mag || curr_mag == 0
            continue;
        end
        %determine radius of object
        radius = round(-alpha*(curr_mag - max_mag) + min_r);
        
        %ints = mag_to_int(curr_mag, band); %calcuate the intesity at pix
        
        %create submask of object (create array of 0's and 1's,
        %where the 1's represnt to location of the objects in the
        %submask). Basically a sqaure of 0's with a circle of 1's
        [rr, cc] = meshgrid(1:2*radius+1);
        Circle = sqrt((rr-radius-1).^2+(cc-radius-1).^2)<=radius;
        
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
                if Circle(i,j) == 1
                    Z(ycurr,xcurr) = NaN;
                end
            end
        end
    end
end

%%write masked image
%newfits = sprintf('%s_masksubtracted.fits', band);
%imagedata = single(Z);
%fitswrite(imagedata, newfits);
end


