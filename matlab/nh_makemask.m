
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_makemask.m
%%  Jun 2016, MZ
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,mu] = nh_makemask(data,paths,params,nsig,use_gaia,new_star_mask,max_mag,save_file,flag_method,errflag_mags)

% load('run_params.mat','params')
if params.err_mags == 1
    % Choose data dir based on err_flag
    datastruct = data.(params.err_str);
else
    datastruct = data;
end

% read in the .fits data for the calibrated image file
imagename = ['regist_',datastruct.header.rawfile];
image_i = fitsread(sprintf('%s%s',paths.imagedir,imagename));

% read in datastruct.data to replace possibly modified version from earlier run
datastruct.data = image_i;

%% here begins the star catalog masking

% load up the corresponding catalog file
if use_gaia == 1
    load(sprintf('%smat_files/field_%d_data.mat',paths.gaiadir,datastruct.header.fieldnum));
    
    % figure out the length of the catalog
    [ncat,~] = size(RA);
elseif use_gaia == 0
    load(sprintf('%sfield_%d_data.mat',paths.catdir,datastruct.header.fieldnum));
    
    % figure out the length of the catalog
    [~,ncat] = size(RA);
end

% figure out the size the mask arrays need to be
xdim = datastruct.astrom.imagew;
ydim = datastruct.astrom.imageh;

% make the empty mask arrays
fullmask = zeros(xdim,ydim);
ghostlessmask = zeros(xdim,ydim);
linemask = zeros(xdim,ydim);
starmask = zeros(xdim,ydim);
statmask = zeros(xdim,ydim);
clipmask = zeros(xdim,ydim);

% the maximum magnitude of star to mask - this should be a value that
% gives back approximately 50% of pixels masked, and as large as possible
if strcmp(flag_method,'new') == 1
    max_mag = max_mag; %old way: 16.5 + 2.5 * log10(sqrt(datastruct.astrometry.id_exptime));
elseif strcmp(flag_method, 'old') == 1
    max_mag = 17.75;
elseif strcmp(flag_method, 'old_corr') == 1
    max_mag = 16.5 + 2.5 * log10(sqrt(datastruct.astrometry.id_exptime));
end

%%%%%%%%%%%%%%%% old stuff
%Set values for radius calculation: r(mag) = -alpha(mag - max_mag)+min_r
%min_r = 1;6
%alpha = 0.75; %needs to be calculated
%beta = 2;
%bright_mag = 3;
%%%%%%%%%%%%%%%

% this will hold the
mymag = zeros(ncat,1);
myflux = zeros(ncat,1);
myxpix = zeros(ncat,1);
myypix = zeros(ncat,1);
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

if new_star_mask == 1 || ~isfield(datastruct,'mask')
    
    % Create star mask - takes very long and can skip if already saved
    % loop over each catalog entry
    for row = 1:ncat
        
        %if using USNOB1 catalog
        if use_gaia == 0
            % first, we have this mish-mash set of magnitudes.  They seem to allow
            % either 0 to indicate non-detection(?) or crazy high values which are
            % below the survey limit.
            % so to deal with this, I have decided to use only t6hose sources within
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
            [ypix, xpix] = radec2pix(RA(row),DEC(row), datastruct.astrom);
            
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
                        linemask(:,round(xpix)) = 1; % this is set to mask column, not row
                        % so far, this isn't used in testing set or for lauer
                    end

                    if thismag < 10 & xpix >= 1 && xpix <= xdim && ypix >= 1 && ypix <= ydim
                        linemask(round(ypix),:) = 1; % this does row
                    end
                    
                    % another little piece of housekeeping; just making sure that we're
                    % keeping track of the magnitudes of stars that have made it this far.
                    mymag(row) = thismag;
                    myxpix(row) = xpix;
                    myypix(row) = ypix;
                    
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
            
            %if using gaia catalog
        elseif use_gaia == 1
            
            %mags don't seem to have any wild values, so we'll use all of
            %them without restriction
            if errflag_mags == 1
                thismag = Gmag(row) + Gmagerr(row)*randn(1); %If including mag error, include up to the full reported Gaia error with Gaussian probability
            elseif errflag_mags == 0
                thismag = Gmag(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value
            end
                
            allmags(row) = thismag;
            
            %find x/y coordinate of the object
            [ypix, xpix] = radec2pix(RA(row),DEC(row), datastruct.astrom);
            
            % prepare to sum flux of object being masked
            fluxsum = 0;
            % check if the object is in the image
            if xpix >= 1-20 && xpix <= xdim+20 && ypix >= 1-20 && ypix <= ydim+20
                numinimg = numinimg + 1;
                
                % require that the magnitude is within sensible bounds
                if thismag < max_mag & ~isnan(thismag)
                    numinbnds = numinbnds + 1;
                    
                    % this is responding to the observation that bright sources leave
                    % charge transfer tracks in these images.  These must be masked.  For
                    % sources brighter than V=7, we mask the row the star falls into.
                    if thismag < 7 & xpix >= 1 && xpix <= xdim && ypix >= 1 && ypix <= ydim
                        linemask(:,round(xpix)) = 1; % this is set to mask column, not row
                        % so far, this isn't used in testing set or for lauer
                    end
                        
                    % 13 mag seems to catch all stars that have artifacts
                    if thismag < 13 & xpix >= 1 && xpix <= xdim && ypix >= 1 && ypix <= ydim
                        linemask(round(ypix),round(xpix):end) = 1; % this does row from only middle of star to right side
                    end
                    
                    % another little piece of housekeeping; just making sure that we're
                    % keeping track of the magnitudes of stars that have made it this far.
                    %                     mymag(row) = thismag;
                    %                     myxpix(row) = xpix;
                    %                     myypix(row) = ypix;
                    
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
                                fluxsum = fluxsum + datastruct.data(ycurr,xcurr)./ datastruct.astrom.exptime;
                            end
                        end
                    end
                    
                    if fluxsum ~= 0
                        mymag(row) = thismag;
                        myxpix(row) = xpix;
                        myypix(row) = ypix;
                        myflux(row) = fluxsum;
                    end
                    
                    %checks for why things are being excluded and how many
                elseif thismag > max_mag
                    checkmag(row) = thismag;
                    nummagmax = nummagmax + 1;
                elseif isnan(thismag)
                    nummagnan = nummagnan + 1;
                end
            end
        end
    end
    
    %     min(mymag)
    myxpix = myxpix(myxpix~=0);
    myypix = myypix(myypix~=0);
    
    if save_file == 1
        %keeping running lists of source mags, coordinates, and fluxes
        mymag = mymag(mymag~=0);
        myflux = myflux(myflux~=0);
        mycalflux = datastruct.cal.sbconv .* myflux;
        %save list of stars to file
        %     star_list = horzcat(myxpix,myypix,mymag);
        star_list = horzcat(mymag,mycalflux);
        
        filename = strcat(num2str(max_mag),'_mask.mat');
        save(filename,'star_list');
        %     filename = strcat(paths.starlistdir,datastruct.header.timestamp,'star_list.mat');
        %     save(filename,'star_list');
    end
    
elseif new_star_mask == 0
    
    % If not repeating star mask creation, read in from data file
    starmask = datastruct.mask.starmask;
    linemask = datastruct.mask.linemask;
    
end
%%%%%%%%%%%%%%%%%%
%% thus ends the star catalog masking
%% here begins the manual masking

% If analyzing old ghost data, don't manually mask (masks out ghosts!)
if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/ghosts/') == 1
    manmask = zeros(size(datastruct.data));
else
    filein = sprintf('%s%s.mat',paths.mandir,datastruct.header.timestamp);
    % If man mask files exists, create the mask, otherwise no man mask
    if numel(dir(filein))
        load(filein);
    else
        manmask = zeros(size(datastruct.data));
    end
end

%%%%%%%%%%%%%%%%%%
%% thus ends the manual masking
%% here begins the optical ghost masking

if strcmp(flag_method,'new') == 1
    % If image has a star > mag 8 in range of causing a ghost, mask the ghost
    if ~isempty(datastruct.ghost.brightmag) %datastruct.ghost.brightmag(1,1) > 0
        [ghostmask,datastruct] = nh_ghostmask(datastruct,paths);
    else
        ghostmask = zeros(256);
    end
else
    ghostmask = zeros(256);
end

%%%%%%%%%%%%%%%%%%
%% thus ends the optical ghost masking
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
    [xp,yp] = meshgrid([1:size(datastruct.data,1)]-256./2);
    rad = sqrt(xp.^2 + yp.^2);
    whpl = rad <= 1.5.*50; % mask 50% futher than necessary based on Cheng et
    %  al description
    statmask(whpl) = 1;
end

%% that's the static mask (right now just clips out edge pixels)
%% here begins the clip masking
%% first, compute the mean and std of the unmasked pixels

if strcmp(flag_method,'new') == 1
    % 3-round clip mask
    for iter=1:3
        if iter == 1
            clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
            clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
            % now find all unmasked pixels > 3 sigma away from the mean
            % This method does +/- 3sig from mean, which could overclip on the negative side since the tail tends to be positive
%             whpl = ((datastruct.data > clipmean + nsig.*clipstd) | (datastruct.data < clipmean - nsig.*clipstd)) &...
%                 ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
            % This method does + 3sig from mean on abs(data) so neg side not clipped as hard as pos
            whpl = (abs(datastruct.data) > clipmean + nsig.*clipstd) &...
                ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
            clipmask(whpl) = 1;
        elseif iter == 2
            clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask));
            clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask));
            % now find all unmasked pixels > 3 sigma away from the mean
            whpl = (abs(datastruct.data) > clipmean + nsig.*clipstd) &...
                ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask;
            clipmask(whpl) = 1;
        elseif iter == 3
            clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask));
            clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask));
            % now find all unmasked pixels > 2 sigma away from the mean
            whpl = (abs(datastruct.data) > clipmean + (nsig-1).*clipstd) &...
                ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask;
            clipmask(whpl) = 1;
        end
    end
    % 2-round clip mask
%     for iter=1:2
%         if iter == 1
%             clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
%             clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
%             % now find all unmasked pixels > 3 sigma away from the mean
%             % This method does +/- 3sig from mean, which could overclip on the negative side since the tail tends to be positive
% %             whpl = ((datastruct.data > clipmean + nsig.*clipstd) | (datastruct.data < clipmean - nsig.*clipstd)) &...
% %                 ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
%             % This method does + 3sig from mean on abs(data) so neg side not clipped as hard as pos
%             whpl = (abs(datastruct.data) > clipmean + nsig.*clipstd) &...
%                 ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
%             clipmask(whpl) = 1;
%         elseif iter > 1
%             clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask));
%             clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask));
%             % now find all unmasked pixels > 3 sigma away from the mean
%             whpl = (abs(datastruct.data) > clipmean + nsig.*clipstd) &...
%                 ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask & ~clipmask;
%             clipmask(whpl) = 1;
%         end
%     end
    
elseif strcmp(flag_method,'old') == 1
    clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
    clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
    
    % old way - seems wrong - instead of masking data >/< mean +/- nsig*std it
    % masks abs(data) > mean + nsig*std, this assumes mean is ~0 and
    % gaussian is symmetric around 0, but mean is positive so mean - nsig*std not
    % the same as mean + nsig*std
    whpl = (abs(datastruct.data) > clipmean + nsig.*clipstd) & ...
        ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
    
    % and set the mask bit
    clipmask(whpl) = 1;
    
elseif strcmp(flag_method,'old_corr') == 1
    
    clipmean = mean(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
    clipstd = std(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask));
    
    % now find all unmasked pixels > 3 sigma away from the mean
    whpl = ((datastruct.data > clipmean + nsig.*clipstd) | (datastruct.data < clipmean - nsig.*clipstd)) &...
        ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
    
    
    %     % Median Option see https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
    %     clipmedian = median(datastruct.data(~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask)); %  Get the median of the data
    %     clipmediandist = abs(clipmedian - datastruct.data); % Get the absolute distance between the data and the median of the data
    %     clipmediandistcompare = clipmediandist/clipmedian; % Scale the dist from the median by the median of the dist from the median
    %     whpl = (clipmediandistcompare > 3.5) &... % Get where the scaled dist is greater than the comparator value
    %         ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
    
    % whpl = (abs(datastruct.data) > 10) & ...
    %     ~starmask & ~linemask & ~statmask & ~manmask & ~ghostmask;
    
    % and set the mask bit
    clipmask(whpl) = 1;
end

% Turn off clip mask
%     clipmask = zeros(size(datastruct.data));

%%%%%%%%%%%%%%%%%
%%  now create the union of all of the masks as a bitwise mask
whpl = starmask == 1;
fullmask(whpl) = 1;
ghostlessmask(whpl) = 1;
whpl = linemask == 1;
fullmask(whpl) = fullmask(whpl) + 2;
ghostlessmask(whpl) = ghostlessmask(whpl) + 2;
whpl = clipmask == 1;
fullmask(whpl) = fullmask(whpl) + 4;
ghostlessmask(whpl) = ghostlessmask(whpl) + 4;
whpl = statmask == 1;
fullmask(whpl) = fullmask(whpl) + 8;
ghostlessmask(whpl) = ghostlessmask(whpl) + 8;
whpl = manmask == 1;
fullmask(whpl) = fullmask(whpl) + 16;
ghostlessmask(whpl) = ghostlessmask(whpl) + 16;
whpl = ghostmask == 1;
fullmask(whpl) = fullmask(whpl) + 32;

% and so that we have it to hand, compute the mean of the fully masked pixels
datmean = mean(datastruct.data(~fullmask)./ datastruct.astrom.exptime); %DN/s
datstd = std(datastruct.data(~fullmask)./ datastruct.astrom.exptime);

% create a mask which is just 1 everywhere there is a bad pixel, no matter
% the reason
onemask = fullmask > 0;
oneghostlessmask = ghostlessmask > 0; % mask without ghostmask, used in ghost_analysis for diff ghost fit 

% % Calculate peak hist of masked image
% %Preallocate for edge values of peak histogram value
% minmaxedge = zeros(2,1);
% %Save only ghost area of image
% im = datastruct.datastruct.*~datastruct.mask.onemask;
% %Find indices where pixel values are > 0
% idx = im>0;
% %Save the bin edges and bin counts for a histogram of the ghost
% [N,edges] = histcounts(im(idx));
% [~,I] = max(N); % Find bin with max counts
% minmaxedge(1,1) = edges(I); % Get the edge values for that bin
% minmaxedge(2,1) = edges(I+1);
% idx2 = (im >= edges(I)) & im < (edges(I+1)); % Get the values in the bin edges (uses histcount bin edge logic)
% [N,edges] = histcounts(im(idx&idx2)); % Redo histcounts with only edges, let auto alg work
% edges = linspace(edges(1),edges(end),length(edges)*10); % Boost edge fidelity
% %     bin_width = mean(diff(edges))
% [N,edges] = histcounts(im(idx&idx2),edges); % Redo histcounts again with higher fidelity bin edges based on auto alg
% %Find the bin with maximum counts
% [M,I] = max(N);
% %Calculate edge values for that bin
% minmaxedge(1,1) = edges(I);
% minmaxedge(2,1) = edges(I+1);
% %Pixel value is average of those bin edges
% most_prob_val = median(minmaxedge)/datastruct.astrom.exptime;
%     datmean2 = mean(datastruct.data(~fullmask&idx)./ datastruct.astrom.exptime)
%     datmedian = median(datastruct.data(~fullmask&idx)./ datastruct.astrom.exptime)

% compute the mask fraction
maskfrac = sum(onemask(:))./256.^2;

% dump all of this information into a structure
mask.mask = fullmask;
mask.starmask = starmask;
mask.linemask = linemask;
mask.clipmask = clipmask;
mask.statmask = statmask;
mask.manmask = manmask;
mask.ghostmask = ghostmask;
mask.oneghostlessmask = oneghostlessmask;
mask.onemask = onemask;
mask.maskfrac = maskfrac;
mask.maxmag = max_mag;

% If error is on, write to error substruct
if params.err_mags == 1
    datastruct.(params.err_str).mask = mask;
    datastruct.(params.err_str).stats.maskmean = datmean;
    datastruct.(params.err_str).stats.maskstd = datstd;
    datastruct.(params.err_str).stats.maskerr = datstd ./ sqrt(256.^2 - sum(onemask(:)));
else
    % and append it to the data structure
    datastruct.mask = mask;
    datastruct.stats.maskmean = datmean;
    datastruct.stats.maskmean_dns_frommask = datmean;
    % datastruct.stats.mostprob = most_prob_val;
    datastruct.stats.maskstd = datstd;
    datastruct.stats.maskerr = datstd ./ sqrt(256.^2 - sum(onemask(:)));
end

%Plot masked data and non-masked data on same plot
% ax1 = subplot(1,2,1);
% imagesc(datastruct.datastruct.*~datastruct.mask.onemask)
% colorbar
% caxis([-5,100])
% caxis('auto')

% ax2 = subplot(1,2,2);
% imagesc(datastruct.data)
% colorbar
% caxis([-5,100])
% caxis(ax1.CLim)

% %Plot masked data with optional mouse movement returns value
h = figure();
clf;
imagesc(datastruct.data.*~datastruct.mask.onemask)
set(h,'visible','off');
% set (gcf, 'WindowButtonMotionFcn', @mouseMove);
a = colorbar;
a.Label.String = 'Intensity [DN]';
pbaspect([1 1 1]);
xlabel('LORRI X Pixels');
ylabel('LORRI Y Pixels');
caxis([-10,10]);
%         title(sprintf('%s',datastruct.header.rawfile));
% grid minor;
% title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
title(sprintf('Field: %d',datastruct.header.fieldnum));
set(gca,'YDir','normal');
ext = '.png';
if not(isfolder([paths.maskdir]))
    mkdir([paths.maskdir])
end
imagename = sprintf('%s%s%s',paths.maskdir,datastruct.header.timestamp,ext);
% imagename = sprintf('%s%s%s',paths.selectmaskdir,datastruct.header.timestamp,ext);
print(h,imagename, '-dpng');

% Plot histogram of masked image using hist fit and overplotting mean +/-
% sigma
maskim = datastruct.data.*~datastruct.mask.onemask;
maskim(maskim==0) = NaN;
bins = 40;

h = figure();
clf;
set(h,'visible','off');
histfitCustom(maskim(:),bins,'normal');
title(['Mean = ',num2str(mean(maskim(:),'omitnan')),' Median = ',num2str(median(maskim(:),'omitnan'))]);
ylimMax = ylim; %get the ylims (min and max)
ylimMax = ylimMax(2); %get the actual max
ylim( [0.5, ylimMax] ); %force 0
set(gca,'YScale','log');
hold on
% Mean +/- sig
%     plot( repmat(nanmean(maskim(:)),5,1),linspace(0.5,ylimMax,5),'linewidth',5,'color','green');
%     plot( repmat(nanmean(maskim(:))+nanstd(maskim(:)),5,1),linspace(0.5,ylimMax,5),'linewidth',5,'color','green');
%     plot( repmat(nanmean(maskim(:))-nanstd(maskim(:)),5,1),linspace(0.5,ylimMax,5),'linewidth',5,'color','green');
% Mean and median
plot( repmat(mean(maskim(:),'omitnan'),5,1),linspace(0.5,ylimMax,5),'linewidth',1,'color','green');
plot( repmat(median(maskim(:),'omitnan'),5,1),linspace(0.5,ylimMax,5),'linewidth',1,'color','magenta');
hold off
legend('','','Mean','Median')
xlabel('Unmasked Pixel Intensity [DN]');
ylabel('N');
ext = '.png';
if not(isfolder([paths.histdir]))
    mkdir([paths.histdir])
end
imagename = sprintf('%s%s%s',paths.histdir,datastruct.header.timestamp,ext);
print(h,imagename, '-dpng');

%Calculate fit to data
pd = fitdist(maskim(:),'normal');
mu = pd.mu;
datmean = mean(datastruct.data(~fullmask));
%Calculate the theoretical PDF from the fit parameters
x_range = linspace(floor(nanmin(nanmin(maskim))),ceil(nanmax(nanmax(maskim))),100);
probability_predicted = pdf(pd,x_range);

%     % Plot histogram with calculated fit over the top
%     figure
%     % histogram(maskim(:),60,'Normalization','probability')
%     g = histogram(maskim(:),bins);
%     hold on
%     h = plot(x_range,probability_predicted*max(g.Values)/max(probability_predicted),'.-','LineWidth',3);
%     set(gca,'YScale','log');
%     ylimMax = ylim; %get the ylims (min and max)
%     ylimMax = ylimMax(2); %get the actual max
%     ylim( [0.5, ylimMax] ); %force 0
%     xlim([-25,25]);

%Plot histfit and calculated fit together
%     h = figure;
%     clf;
%     set(h,'visible','off');
%     % histfit(maskim(:),bins,'normal');
%     g = histfitCustom(maskim(:),bins,'normal');
%
%     hold on
% %     plot(x_range,probability_predicted*max(g(1).YData)/max(probability_predicted),'.-','LineWidth',3);
%     set(gca,'YScale','log');
%     title(['\mu fit = ',num2str(pd.mu)]); %sprintf('%s',datastruct.header.rawfile)
%     ylim( [0.5, 100000] ); %force 0
% %     xlim([-50,50]);
%     ext = '.png';
%     imagename = sprintf('%s%s%s',paths.histdir,datastruct.header.timestamp,ext);
%     print(h,imagename, '-dpng');


%inspired by https://www.mathworks.com/matlabcentral/answers/313862-how-to-plot-a-normal-distribution-graph-to-fit-a-bar-graph
mu = nanmean(maskim(:));
%     sd = nanstd(maskim(:));
%     fun_dist_norm = @(mu,sd,x) exp(-(x-mu).^2 ./ (2*sd^2)) /(sd*sqrt(2*pi)); % Function for Standard Normal Distribution
%
%     [counts,edges] = histcounts(maskim(:), bins);  % Histogram
%
%     centers = edges(1:length(edges)-1) + mean(diff(edges))/2; % Calculate centers
%
%     centers_xVect = linspace(min(centers), max(centers),500); % x-vector for plotting
%
%     norm_dist = fun_dist_norm(mu,sd,centers_xVect); % Calculate Standard Normal Distribution
%
%     h = figure;
%     bar(centers, counts,1); % Plot Histogram
%     hold on
%     plot(centers_xVect, norm_dist*max(counts)/max(norm_dist), '.-', 'LineWidth',2); % Plot Scaled Standard Normal Distribution
%     hold off
%     set(gca,'YScale','log');
%     ylimMax = ylim; %get the ylims (min and max)
%     ylimMax = ylimMax(2); %get the actual max
%     ylim( [0.5, ylimMax] );
%     xlim([-25,25]);
%     grid
%see
%https://www.mathworks.com/matlabcentral/answers/175958-function-pdf-doesn-t-return-pdf-values
%as to why norm_dist (or histfit for that matter) don't return perfectly
%scaled results


if params.err_mags == 1
    % Choose data dir based on err_flag
    data.(params.err_str) = datastruct;
else
    data = datastruct;
end

end
