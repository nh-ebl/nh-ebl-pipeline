% This function reads in USNOB1 catalog files for each image and searches
% for stars with mag<8 between 0.145 and 0.37 degrees from boresight
% Star magnitudes and x/y coordinates are recorded
% Set up to be run from nh_pipeline

% Symons, May 2019
% Thayer, June 2019
%added the RA and DEC variables and recorded their values

function [data, ghostcount, realghostcount] = nh_findghoststar(data, paths, use_gaia)

% load up the corresponding catalog file
if use_gaia == 1
    load(sprintf('%smat_files/field_%d_data.mat',paths.gaiadir,data.header.fieldnum));

    % figure out the length of the catalog
    [ncat,~] = size(RA);
elseif use_gaia == 0
    load(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum));

    % figure out the length of the catalog
    [~,ncat] = size(RA);
end

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
brightra=zeros(ncat,1);
brightdec=zeros(ncat,1);
raa=zeros(ncat,1);
decc=zeros(ncat,1);
boxxpix = zeros(ncat,1);
boxypix = zeros(ncat,1);
boxmag = zeros(ncat,1);

% loop over each catalog entry;
for row = 1:ncat
    
    %if using USNOB1 catalog
    if use_gaia == 0
    
        % retrieve magnitudes and star/galaxy classification
        mags = [B1mag(row),B2mag(row),R1mag(row),R2mag(row),I2mag(row)];
        sg = [B1sg(row),B2sg(row),R1sg(row),R2sg(row),I2sg(row)];
        raa=RA(row);
        decc=DEC(row);
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
        % if all star/galaxy classfifications are nan (signifying bright star),
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
    
    elseif use_gaia == 1
        %mags don't seem to have any wild values, so we'll use all of
        %them without restriction
        thismag = Gmag(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value
        raa=RA(row);
        decc=DEC(row);

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
            % save star to list of stars in ghost box
            boxxpix(row) = xpix;
            boxypix(row) = ypix;
            boxmag(row) = thismag;
            
            % require that the magnitude is within sensible bounds
            if thismag < 8 && ~isnan(thismag)
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
            end
        end
    end
end

% remove zero values
quadxpix = quadxpix(quadxpix~=0);
quadypix = quadypix(quadypix~=0);
quadmag = quadmag(quadmag~=0);
brightra= raa(raa~=0);
brightdec=decc(decc~=0);
boxxpix = boxxpix(boxxpix~=0);
boxypix = boxypix(boxypix~=0);
boxmag = boxmag(boxmag~=0);

if numinbnds > 0
    ghostcount = 1;
else
    ghostcount = 0;
end

% check for logged physical ghost - requires ghost_analysis to have already
% been run and saved to data - if no data.ghost yet, can't do this
if strcmp(data.ghost.ghostx , '') == 1 || data.ghost.ghostx == 0
    realghostcount = 0;
else
    realghostcount = 1;
end
% set count to 0 if no logged ghost data
% realghostcount=0;

% disp(min(Gmag))

% save bright star info to data
data.ghost.brightmag = quadmag;
data.ghost.brightxpix = quadxpix;
data.ghost.brightypix = quadypix; 
data.ghost.brightra = raa;
data.ghost.brightdec = decc;
data.ghost.boxxpix = boxxpix;
data.ghost.boxypix = boxypix;
data.ghost.boxmag = boxmag;
end