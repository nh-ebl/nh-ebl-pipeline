% This function reads in USNOB1 catalog files for each image and searches
% for stars with mag<8 between 0.145 and 0.37 degrees from boresight
% Star magnitudes and x/y coordinates are recorded
% Set up to be run from nh_pipeline

% Symons, May 2019
% Thayer, June 2019
%added the RA and DEC variables and recorded their values

function [data, ghostcount, realghostcount] = nh_findghoststar(data, paths, params, use_gaia, errflag_mags, mc)

% If doing err gals MC, load appropriate catalog, if not, load regular catalog
if params.err_gals == 0
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
elseif params.err_gals == 1
    % load up the corresponding catalog file
    if use_gaia == 1
        load(sprintf('%smat_files/field_%d_mc/%i', paths.gaiadir,data.header.fieldnum,mc));
        
        % figure out the length of the catalog
        [ncat,~] = size(RA);
    end
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

if ~exist('xpix','var') || ~exist('ypix','var')
    %only run if need to make xpix/ypix (should not run since makemask made them)
    ypix = zeros(ncat,1);
    xpix = zeros(ncat,1);
    % convert radec2pix parallel now since will be required (later loop
    % can't be paralleled w/o work)
    parpoolobj = gcp('nocreate'); % check for thread pool
    if isempty(parpoolobj)
        maxNumCompThreads(7); % Declare number of threads to use
        parpool('threads');
    else
        if ~isa(parpoolobj,'parallel.ThreadPool')
            delete(parpoolobj); %want threads here
            maxNumCompThreads(7); % Declare number of threads to use
            parpool('threads');
        end
    end
    RApar = RA;
    DECpar = DEC; %parfor was mad about these not ever being declared officially
    parfor row = 1:ncat
        %parallel for speed, later loop needs work to parallize
        [ypix(row), xpix(row)] = radec2pix(RApar(row),DECpar(row), data.astrom);
    end
    if params.err_gals == 0
        % save the calc'd ypix/xpix the corresponding catalog file
        if use_gaia == 1
            save(sprintf('%smat_files/field_%d_data.mat',paths.gaiadir,data.header.fieldnum),'xpix','ypix','-append');
        elseif use_gaia == 0
            save(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum),'xpix','ypix','-append');
        end
    elseif params.err_gals == 1
        % save the calc'd ypix/xpix the corresponding catalog file
        if use_gaia == 1
            save(sprintf('%smat_files/field_%d_mc/%i', paths.gaiadir,data.header.fieldnum,mc),'xpix','ypix','-append');
        end
    end
end

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
        % [ypix(row), xpix(row)] = radec2pix(RA(row),DEC(row), data.astrom); %par now
        
        % check if the object is in the target zone
        check = 0;
        % quad 1
        if xpix(row) >= 1-199 && xpix(row) <= 1 && ypix(row) >= 1-199 && ypix(row) <= ydim+199
            quad1 = quad1 + 1;
            check = 1;
        end
        % quad 2
        if xpix(row) >= xdim && xpix(row) <= xdim + 199 && ypix(row) >= 1-199 && ypix(row) <= ydim+199
            quad2 = quad2 + 1;
            check = 1;
        end
        % quad 3
        if xpix(row) >= 1 && xpix(row) <= xdim && ypix(row) >= ydim && ypix(row) <= ydim+199
            quad3 = quad3 + 1;
            check = 1;
        end
        % quad 4
        if xpix(row) >= 1 && xpix(row) <= xdim && ypix(row) >= 1-199 && ypix(row) <= 1
            quad4 = quad4 + 1;
            check = 1;
        end
        
        if check == 1
            % require that the magnitude is within sensible bounds
            if thismag < 8 && ~isnan(thismag) && thissg > 0
                numinbnds = numinbnds + 1;
                
                % save x/y location of star and magnitude
                quadxpix(row) = xpix(row);
                quadypix(row) = ypix(row);
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
        if errflag_mags == 1
            thismag = Gmag(row) + Gmagerr(row)*randn(1); %If including mag error, include up to the full reported Gaia error with Gaussian probability
        elseif errflag_mags == 0
            thismag = Gmag(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value
        end
        raa=RA(row);
        decc=DEC(row);
        
        %find x/y coordinate of the object
        % [ypix(row), xpix(row)] = radec2pix(RA(row),DEC(row), data.astrom); %par now
        
        % check if the object is in the target zone
        check = 0;
        % quad 1
        if xpix(row) >= 1-199 && xpix(row) <= 1 && ypix(row) >= 1-199 && ypix(row) <= ydim+199
            quad1 = quad1 + 1;
            check = 1;
        end
        % quad 2
        if xpix(row) >= xdim && xpix(row) <= xdim + 199 && ypix(row) >= 1-199 && ypix(row) <= ydim+199
            quad2 = quad2 + 1;
            check = 1;
        end
        % quad 3
        if xpix(row) >= 1 && xpix(row) <= xdim && ypix(row) >= ydim && ypix(row) <= ydim+199
            quad3 = quad3 + 1;
            check = 1;
        end
        % quad 4
        if xpix(row) >= 1 && xpix(row) <= xdim && ypix(row) >= 1-199 && ypix(row) <= 1
            quad4 = quad4 + 1;
            check = 1;
        end
        
        if check == 1
            % save star to list of stars in ghost box
            boxxpix(row) = xpix(row);
            boxypix(row) = ypix(row);
            boxmag(row) = thismag;
            
            % require that the magnitude is within sensible bounds
            if thismag < 8 && ~isnan(thismag)
                numinbnds = numinbnds + 1;
                
                % save x/y location of star and magnitude
                quadxpix(row) = xpix(row);
                quadypix(row) = ypix(row);
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
if isfield(data,'ghost')
    if ~isfield(data.ghost,'ghostx')
        realghostcount = 0;
    elseif strcmp(data.ghost.ghostx , '') == 1 || data.ghost.ghostx == 0
        realghostcount = 0;
    else
        realghostcount = 1;
    end
else
    realghostcount=0;
end
% set count to 0 if no logged ghost data
% realghostcount=0;

% disp(min(Gmag))

% If error is on, write to error substruct
% load('run_params.mat','params')
if params.err_mags == 1 || params.err_gals == 1
    % save bright star info to data
    data.(params.err_str).ghost.brightmag = quadmag;
    data.(params.err_str).ghost.brightxpix = quadxpix;
    data.(params.err_str).ghost.brightypix = quadypix;
    data.(params.err_str).ghost.brightra = raa;
    data.(params.err_str).ghost.brightdec = decc;
    data.(params.err_str).ghost.boxxpix = boxxpix;
    data.(params.err_str).ghost.boxypix = boxypix;
    data.(params.err_str).ghost.boxmag = boxmag;
else
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

end