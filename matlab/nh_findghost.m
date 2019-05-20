% This function reads in USNOB1 catalog files for each image and searches
% for stars with mag<8 between 0.145 and 0.37 degrees from boresight
% Star magnitudes and x/y coordinates are recorded
% Set up to be run from nh_pipeline

% Symons, May 2019

function data = nh_findghost(data)

% get the paths for this file
paths = get_paths();

% load up the corresponding catalog file
load(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum));

% figure out the length of the catalog
[~,ncat] = size(RA);

% figure out the size the mask arrays need to be
xdim = data.astrom.imagew;
ydim = data.astrom.imageh;

% preallocate
numinimg = 0;
numinbnds = 0;
checkmag = zeros(ncat,1);
checksg = zeros(ncat,1);
rangeexclude = 0;
nummagmax = 0;
numsgover = 0;
nummagnan = 0;
nansgcnt = 0;
quad1 = 0;
quad2 = 0;
quad3 = 0;
quad4 = 0;
quadxpix = zeros(ncat,1);
quadypix = zeros(ncat,1);
quadmag = zeros(ncat,1);

% loop over each catalog entry;
for row = 1:ncat
    
    % retrieve magnitudes and star/galaxy classification
    mags = [B1mag(row),B2mag(row),R1mag(row),R2mag(row),I2mag(row)];
    sg = [B1sg(row),B2sg(row),R1sg(row),R2sg(row),I2sg(row)];
    
    % using all magnitudes, calculate average magnitude
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
    thissg = nanmean(sg);
    % if all star/galaxy classifications are nan (signifying bright star),
    % assign sg of 1 so bright stars not excluded
    if isnan(B1sg(row)) && isnan(B2sg(row)) && isnan(R1sg(row)) && isnan(R2sg(row)) && isnan(I2sg(row))
        nansgcnt = nansgcnt + 1;
        thissg = 1;
    end
    
    %find x/y coordinate of the object
    [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
    
    % check if the object is in the target zone
    check = 0;
    % quad 1
    if xpix >= 1-199 && xpix <= 1 && ypix >= 1-199 && ypix <= ydim+199
        quad1 = quad1 + 1;
        check = 1;
    end
    % quad 2
    if xpix >= xdim && xpix <= xdim + 199 && ypix >= 1-199 && ypix <= ydim+199
        quad2 = quad2 + 1;
        check = 1;
    end
    % quad 3
    if xpix >= 1 && xpix <= xdim && ypix >= ydim && ypix <= ydim+199
        quad3 = quad3 + 1;
        check = 1;
    end
    % quad 4
    if xpix >= 1 && xpix <= xdim && ypix >= 1-199 && ypix <= 1
        quad4 = quad4 + 1;
        check = 1;
    end
        
    if check == 1
        % require that the magnitude is within sensible bounds
        if thismag < 8 && ~isnan(thismag) && thissg > 0
            numinbnds = numinbnds + 1;
            
            % save x/y location of star and magnitude
            quadxpix(row) = xpix;
            quadypix(row) = ypix;
            quadmag(row) = thismag;
        
        % checks to see how many stars excluded for which criteria
        elseif thismag >= 8
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

% remove zero values
quadxpix = quadxpix(quadxpix~=0);
quadypix = quadypix(quadypix~=0);
quadmag = quadmag(quadmag~=0);

% save bright star info to data
data.ghost.brightmag = quadmag;
data.ghost.brightxpix = quadxpix;
data.ghost.brightypix = quadypix; 

end