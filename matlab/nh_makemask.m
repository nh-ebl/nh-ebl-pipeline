%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_makemask.m
%%  Jun 2016, MZ
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = nh_makemask(data)

% get the paths for this file
paths = get_paths();

% load up the corresponding catalog file
load(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum));

% figure out the length of the catalog
[~,ncat] = size(RA);

% figure out the size the mask arrays need to be
xdim = data.astrom.imagew;
ydim = data.astrom.imageh;

% make the empty mask arrays
fullmask = zeros(xdim,ydim);
linemask = zeros(xdim,ydim);
starmask = zeros(xdim,ydim);
statmask = zeros(xdim,ydim);
clipmask = zeros(xdim,ydim);

% the maximum magnitude of star to mask - this should be a value that
% gives back approximately 50% of pixels masked, and as large as possible
max_mag = 16.5 + 2.5 * log10(sqrt(data.astrometry.id_exptime));

%%%%%%%%%%%%%%%% old stuff
%Set values for radius calculation: r(mag) = -alpha(mag - max_mag)+min_r
%min_r = 1;
%alpha = 0.75; %needs to be calculated
%beta = 2;
%bright_mag = 3;
%%%%%%%%%%%%%%%

% this will hold the
mymag = zeros(ncat,1);
allmags = zeros(ncat,1);
numinimg = 0;
numinbnds = 0;
checkmag = zeros(ncat,1);
checksg = zeros(ncat,1);
rangeexclude = 0;
nummagmax = 0;
numsgover = 0;
nummagnan = 0;
nansgcnt = 0;

% loop over each catalog entry;
for row = 1:ncat
    
    % first, we have this mish-mash set of magnitudes.  They seem to allow
    % either 0 to indicate non-detection(?) or crazy high values which are
    % below the survey limit.
    % so to deal with this, I have decided to use only those sources within
    % a reasonable limit (1 and 21 mags, the quoted V-band depth of the
    % survey) and then to take the median
    % magnitude between whichever members of the four colors are sane.
    
    mags = [B1mag(row),B2mag(row),R1mag(row),R2mag(row),I2mag(row)];
    sg = [B1sg(row),B2sg(row),R1sg(row),R2sg(row),I2sg(row)];
    
    lambda_mag = [425,462.5,645,650,810];
    whpl = (mags < 20) & (mags > 1);
    if sum(whpl) > 1
        thismag = ...
            nh_synthetic_photometry(lambda_mag(whpl),mags(whpl),'LORRI');%+...
        %randn(1) .* 0.25;
    else
        thismag = NaN;
        rangeexclude = rangeexclude+1;
    end
    allmags(row) = thismag;
    thissg = nanmean(sg);
    if isnan(B1sg(row)) && isnan(B2sg(row)) && isnan(R1sg(row)) && isnan(R2sg(row)) && isnan(I2sg(row))
        nansgcnt = nansgcnt + 1;
        thissg = 1;
    end
    
    %find x/y coordinate of the object
    [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
    
    % check if the object is in the image
    if xpix >= 1-20 && xpix <= xdim+20 && ypix >= 1-20 && ypix <= ydim+20
        numinimg = numinimg + 1;
        
        % require that the magnitude is within sensible bounds
        if thismag < max_mag & ~isnan(thismag) & thissg > 0
            numinbnds = numinbnds + 1;
            
            % this is responding to the observation that bright sources leave
            % charge transfer tracks in these images.  These must be masked.  For
            % sources brighter than V=7, we mask the row the star falls into.
            if thismag < 7 & xpix >= 1 && xpix <= xdim && ypix >= 1 && ypix <= ydim
                linemask(:,round(xpix)) = 1;
            end
            
            % another little piece of housekeeping; just making sure that we're
            % keeping track of the magnitudes of stars that have made it this far.
            mymag(row) = thismag;
            
            %determine radius of object
            %             radius = round(-alpha.*(thismag - max_mag) + min_r) + -beta.*(thismag-bright_mag));
            %             radius = 6000.*exp(-thismag)+2.5;
            %             radius = -1.*(thismag - max_mag) + 2;
            radius = 2.5*(max_mag./thismag).^(2); %+ 2;
            if radius > 120
                dbstop
            end
            %ints = mag_to_int(curr_mag, band); %calcuate the intesity at pix
            
            %create submask of object (create array of 0's and 1's,
            %where the 1's represnt to location of the objects in the
            %submask). Basically a sqaure of 0's with a circle of 1's
            [rr, cc] = meshgrid(1:2*radius+1);
            Circle = sqrt((rr-radius-1).^2+(cc-radius-1).^2)<=radius;
            
            %combined the submask (C) and mask (Z) where xpix, ypix is
            %the center of the object
            for i = 1:(2*radius+1)
                xcurr = round(xpix-radius-1+i);
                if xcurr < 1 || xcurr > xdim;
                    continue;
                end
                for j = 1:(2*radius+1)
                    ycurr = round(ypix-radius-1+j);
                    if ycurr < 1 || ycurr > ydim;
                        continue;
                    end
                    if Circle(i,j) == 1
                        starmask(ycurr,xcurr) = 1;
                    end
                end
            end
        elseif thismag > max_mag
            checkmag(row) = thismag;
            nummagmax = nummagmax + 1;
        elseif isnan(thismag)
            nummagnan = nummagnan + 1;
        elseif thissg <= 0
            checksg(row) = thissg;
            numsgover = numsgover + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%
%% thus ends the star catalog masking
%% here begins the manual masking

filein = sprintf('%s%s.mat',paths.mandir,data.header.timestamp);

if numel(dir(filein))
    load(filein);
else
    manmask = zeros(size(data.data));
end

%%%%%%%%%%%%%%%%%%
%% thus ends the manual masking
%% here begins the static masking

edgelength = 5;
statmask(1:edgelength,:) = 1;
statmask(256-edgelength+1:256,:) = 1;
statmask(:,1:edgelength) = 1;
statmask(:,256-edgelength+1:256) = 1;
%statmask(128-16:128+16,128-16:128+16) = 1;

if 0
    % this bit is in case we want to find out what happens if we mask hard for
    % optical ghosts
    [xp,yp] = meshgrid([1:size(data.data,1)]-256./2);
    rad = sqrt(xp.^2 + yp.^2);
    whpl = rad <= 1.5.*50; % mask 50% futher than necessary based on Cheng et
    %  al description
    statmask(whpl) = 1;
end

%% that's the static mask (right now just clips out edge pixels)
%% here begins the clip masking
%% first, compute the mean and std of the unmasked pixels
datmean = mean(data.data(~starmask & ~linemask & ~statmask & ~manmask));
datstd = std(data.data(~starmask & ~linemask & ~statmask & ~manmask));

% now find all unmasked pixels > 3 sigma away from the mean
whpl = (abs(data.data) > datmean + 3.*datstd) & ...
    ~starmask & ~linemask & ~statmask & ~manmask;
% and set the mask bit
clipmask(whpl) = 1;

%%%%%%%%%%%%%%%%%
%%  now create the union of all of the masks as a bitwise mask
whpl = starmask == 1;
fullmask(whpl) = 1;
whpl = linemask == 1;
fullmask(whpl) = fullmask(whpl) + 2;
whpl = clipmask == 1;
fullmask(whpl) = fullmask(whpl) + 4;
whpl = statmask == 1;
fullmask(whpl) = fullmask(whpl) + 8;
whpl = manmask == 1;
fullmask(whpl) = fullmask(whpl) + 16;

% and so that we have it to hand, compute the mean of the fully masked pixels
datmean = mean(data.data(~fullmask)./ data.astrom.exptime);
datstd = std(data.data(~fullmask)./ data.astrom.exptime);

% create a mask which is just 1 everywhere there is a bad pixel, no matter
% the reason
onemask = fullmask > 0;

% compute the mask fraction
maskfrac = sum(onemask(:))./256.^2;

% dump all of this information into a structure
mask.mask = fullmask;
mask.starmask = starmask;
mask.linemask = linemask;
mask.clipmask = clipmask;
mask.statmask = statmask;
mask.manmask = manmask;
mask.onemask = onemask;
mask.maskfrac = maskfrac;
mask.maxmag = max_mag;

% and append it to the data structure
data.mask = mask;
data.stats.maskmean = datmean;
data.stats.maskstd = datstd;
data.stats.maskerr = datstd ./ sqrt(256.^2 - sum(onemask(:)));

% ax1 = subplot(1,2,1);
% imagesc(data.data.*~data.mask.onemask)
% colorbar
% caxis([-5,100])
% caxis('auto')
% 
% ax2 = subplot(1,2,2);
% imagesc(data.data)
% colorbar
% caxis([-5,100])
% caxis(ax1.CLim)

h = figure(1);
clf;
imagesc(data.data.*~mask.onemask)
set(h,'visible','off');
% set (gcf, 'WindowButtonMotionFcn', @mouseMove);
colorbar; 
caxis([-10,10]);
title(sprintf('%s',data.header.rawfile));
ext = '.png';
imagename = sprintf('%s%s%s',paths.maskdir,data.header.timestamp,ext);
print(h,imagename, '-dpng');

end
