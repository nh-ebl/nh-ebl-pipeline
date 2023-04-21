%make simulated lorri images for the purpose of estimating flux
%contributions from masking and integrated starlight

clc;
clear variables;
% close all;

FLG_cumsum = 0; %set to 1 to use cumsum instead of sum for the plot

dataSets = struct;
dataSets.sets = {'old','new','lauer'};
% dataSets.sets = {'newest'}; %just new for now
for i=1:length(dataSets.sets)
    dataSets.(dataSets.sets{i}) = struct; %make sub-structs
end
%Fill old struct if old is activated
if any(strcmp(dataSets.sets,'old'))
    dataSets.old.paths = get_paths_old();
    dataSets.old.goodfields = [3,5,6,7];
end
%Fill new struct
if any(strcmp(dataSets.sets,'new'))
    dataSets.new.paths = get_paths_new();
    dataSets.new.goodfields = [5,6,7,8];
end
%Fill newest struct
if any(strcmp(dataSets.sets,'newest'))
    dataSets.newest.paths = get_paths_newest();
    dataSets.newest.goodfields = flip([2,4,6,7]);
end
%Fill lauer struct
if any(strcmp(dataSets.sets,'lauer'))
    dataSets.lauer.paths = get_paths_lauer();
    dataSets.lauer.goodfields = [1,2,3,4,5,6,7];
end

%set pipeline settings
use_gaia = 1;
max_mag = 30;
tri_type = 'gaia'; % Select 'ubvri' for interpolated LORRI mags or 'gaia' for G mags

magBinsStepSize = 1; %sets the step size for the mag bins

%make sure parpool is the right type
parpoolobj = gcp('nocreate'); % check for thread pool, which can't use the load call
if isa(parpoolobj,'parallel.ThreadPool')
    delete(parpoolobj); %threads can't use load and will error out
end

PLOT_color = {'#1F77B4','#FF7F0E','#7256C1','#2CA02C','#D81B60','#656565',...
              '#0099A9','#3E0402','#480094','#FF4545','#004D40','#9A1A00',...
              '#224EB3','#249A7D','#FF45CA','#38FF32','#AF5D03','#00004E','#5F0038'};
uniqueField_tot = 0; %prep cntr
for iset=1:length(dataSets.sets)
    uniqueField_tot = uniqueField_tot + length(dataSets.(dataSets.sets{iset}).goodfields); %accumulate
end
uniqueField_id = 1:1:uniqueField_tot; %get unique field IDs


FLG_firstRun = true; %first run flag for plotting
legend_handle = []; %holds the legend handles as they happen
legend_text = {}; %holds the legend handles as they happen
colorCntr = 1; %count the colors

%% from catalog, retrieve lists of stars with mags and coords
for iset=1:length(dataSets.sets)
    disp(['On ',dataSets.sets{iset},' set, figuring out field numbers for each file.'])
    paths = dataSets.(dataSets.sets{iset}).paths; %set the current paths to be correct
    goodfields = dataSets.(dataSets.sets{iset}).goodfields; %set thec urrent goodfields
    set_name = dataSets.sets{iset}; %get the current set name ('old' or 'new' or ...)
    
    datafiles = dir(sprintf('%s*.mat',paths.datadir)); % get files in the field set
    nfiles = size(datafiles,1);
    
    inputMapHolder = cell([length(goodfields),1]); %preallocate inputMap holder - holds a struct for each field
    
    nfiles_fieldnum = zeros(nfiles,1); %preallocate field num holder
    parfor ifile=1:nfiles %faster to read it in once at the beginning parallel style
        data_par = load(sprintf('%s%s',paths.datadir,datafiles(ifile).name)); %load in 1
        nfiles_fieldnum(ifile) = data_par.data.header.fieldnum; %record fieldnum
    end

    for ifield=1:length(goodfields)
        disp(['On ',set_name,' field ',num2str(goodfields(ifield)),' (',num2str(ifield),'/',num2str(length(goodfields)),').']);
        
        ifile_array = find(nfiles_fieldnum == goodfields(ifield)); %find match indexes
        if( isempty(ifile_array) )
            disp(['ERROR: no field match (',num2str(goodfields(ifield)),') for ',set_name,...
                '. quitting fix it, image_sim_field.m, printing fieldnums for all files:'])
            disp(nfiles_fieldnum)
            quit cancel %quit code quick
        else
            headerCntr = length(ifile_array); %cnts through headers to find good file
            headerBad = true; %flag to keep searching for good file (no bad header)
            while( headerBad == true )
                ifile = ifile_array(headerCntr); %choose first one
                load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
                if isfield(data.header,'bad')
                    if( data.header.bad == 0 )
                        headerBad = false; %not bad, good to go
                    else
                        headerCntr = headerCntr - 1; %try next file
                        disp(['Bad file detected in ',set_name,' set/field ',num2str(goodfields(ifield)),' file ',sprintf('%s%s',paths.datadir,datafiles(ifile).name),' (file ',num2str(ifile),' of ',num2str(nfiles),')'])
                    end
                else
                    headerBad = false; %assume OK b/c header.bad doesn't exist
                end
            end
            disp(['Using ',set_name,' set/field ',num2str(goodfields(ifield)),' file ',sprintf('%s%s',paths.datadir,datafiles(ifile).name),' (file ',num2str(ifile),' of ',num2str(nfiles),')']);
        end

        % figure out the size the mask arrays need to be
        xdim = data.astrom.imagew;
        ydim = data.astrom.imageh;
        
        %---declare inputMap parameters---
        inputMap_upSize = [xdim*10,ydim*10]; %declare virtual size to make the grid on
        inputMap_size = [xdim,ydim]; %declare real size to use
        factor = inputMap_upSize(1)/inputMap_size(1); %get factor diff
        
        %---load in PSF data---
        psfPlaceHolder = fspecial('gaussian',50,0.5); %make a temporary 2D Gaussian to act as a PSF [sums to 1 already]
        psf = data.psf.modelpsf; %!!!put real psf here
        psfSize = size(psf); %get psf size
        psfSizeD2 = [floor(psfSize(1)/2),ceil(psfSize(1)/2);floor(psfSize(2)/2),ceil(psfSize(2)/2)]; %get the floor and ceil of PSF size/2 - saves doing math later
        
        %---!!!put real data here!!!---
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
        
        %preallocate lists of star mag, flux, coords
        mymag = zeros(ncat,1);
        myflux = zeros(ncat,1);
        myxpix = zeros(ncat,1);
        myypix = zeros(ncat,1);
        myradius = zeros(ncat,1);
        
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
%                     rangeexclude = rangeexclude+1;
                end
                allmags(row) = thismag;
                thissg = nanmean(sg);
                if isnan(B1sg(row)) && isnan(B2sg(row)) && isnan(R1sg(row)) && isnan(R2sg(row)) && isnan(I2sg(row))
                    nansgcnt = nansgcnt + 1;
                    thissg = 1;
                end
                
                % prepare to sum flux of object being masked
                fluxsum = 0;
                %find x/y coordinate of the object
                [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
                
                % check if the object is in the image
                if xpix >= 1-20 && xpix <= xdim+20 && ypix >= 1-20 && ypix <= ydim+20
%                     numinimg = numinimg + 1;
                    
                    % require that the magnitude is within sensible bounds
                    if thismag < max_mag & ~isnan(thismag) & thissg > 0
%                         numinbnds = numinbnds + 1;
                        
                        % another little piece of housekeeping; just making sure that we're
                        % keeping track of the magnitudes of stars that have made it this far.
                        %                     mymag(row) = thismag;
                        %                     myxpix(row) = xpix;
                        %                     myypix(row) = ypix;
                        
                        %determine radius of object
                        radius = 2.5*(max_mag./thismag).^(2); %+ 2;
                        if radius > 120
                            dbstop
                        end
                        
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
                                    fluxsum = fluxsum + data.data(ycurr,xcurr)./ data.astrom.exptime;
                                end
                            end
                        end
                        
                        if fluxsum ~= 0
                            mymag(row) = thismag;
                            myxpix(row) = xpix;
                            myypix(row) = ypix;
                            myflux(row) = fluxsum;
                            myradius(row) = radius; %record radius
                        end
                        
                    elseif thismag > max_mag
                        checkmag(row) = thismag;
                        nummagmax = nummagmax + 1;
                    elseif isnan(thismag)
%                         nummagnan = nummagnan + 1;
                    elseif thissg <= 0
                        checksg(row) = thissg;
                        numsgover = numsgover + 1;
                    end
                end
                
                %if using gaia catalog
            elseif use_gaia == 1
                
                %mags don't seem to have any wild values, so we'll use all of
                %them without restriction
                thismag = Gmag(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value
                
                allmags(row) = thismag;
                
                %find x/y coordinate of the object
                [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
                
                % prepare to sum flux of object being masked
                fluxsum = 0;
                % check if the object is in the image
                if xpix >= 1-20 && xpix <= xdim+20 && ypix >= 1-20 && ypix <= ydim+20
                    
                    % require that the magnitude is within sensible bounds
                    if thismag < max_mag & ~isnan(thismag)
                        
                        
                        % another little piece of housekeeping; just making sure that we're
                        % keeping track of the magnitudes of stars that have made it this far.
                        %                     mymag(row) = thismag;
                        %                     myxpix(row) = xpix;
                        %                     myypix(row) = ypix;
                        
                        %determine radius of object
                        
                        radius = 2.5*(max_mag./thismag).^(2); %+ 2;
                        if radius > 120
                            dbstop
                        end
                        
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
                                    fluxsum = fluxsum + data.data(ycurr,xcurr)./ data.astrom.exptime;
                                end
                            end
                        end
                        
                        if fluxsum ~= 0
                            mymag(row) = thismag;
                            myxpix(row) = xpix;
                            myypix(row) = ypix;
                            myflux(row) = fluxsum;
                            myradius(row) = radius; %record radius
                        end
                        
                        %checks for why things are being excluded and how many
                    elseif thismag > max_mag
                        checkmag(row) = thismag;
                    elseif isnan(thismag)
                        nummagnan = nummagnan + 1;
                    end
                end
            end
        end
        
        myxpix = myxpix(myxpix~=0);
        myypix = myypix(myypix~=0);
        mymag = mymag(mymag~=0);
        myflux = myflux(myflux~=0);
        myradius = myradius(myradius~=0);
        myflux = (data.cal.vzero.*10.^(-mymag/2.5)).*data.const.Jy.*data.const.nW.*data.cal.nu/((data.cal.pixsize^2).*data.astrom.imagew.*data.astrom.imageh); %convert mymag to myflux instead of using the real-image estimated flux
        %total source surface brightness - multiply by psf to get source with
        %appropriate surface brightness
        
        %% Retrieve ISL sources
        fieldnum = data.header.fieldnum;
        
        % step 1: compute the area of a LORRI image and pull out that many stars
        % from the list at random
        surveyarea = data.astrom.imagew.*data.astrom.imageh.*...
            data.cal.pixsize_arcsec.^2 ./ 3600.^2;
        
%         if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/ghosts/') == 1
%             ghost = 1;
%             old = 0;
%             new = 0;
%         elseif strcmp(paths.datadir,'/data/symons/NH_old_data/mat/good/') == 1
%             old = 1;
%             ghost = 0;
%             new = 0;
%         elseif strcmp(paths.datadir,'/data/symons/nh_data/mat/') == 1
%             new = 1;
%             old = 0;
%             ghost = 0;
%         end
        % load the appropriate trilegal catalog files
        if strcmp(set_name,'old')
            load(sprintf('/home/symons/isl_trilegal/%s/isltrilegal_%02d.mat',tri_type,fieldnum));
        elseif strcmp(set_name,'new')
            load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/trilegal/%s_data/%s/isltrilegal_%02d.mat',set_name,tri_type,fieldnum));
        elseif strcmp(set_name,'newest')
            load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/trilegal/%s_data/%s/isltrilegal_%02d.mat',set_name,tri_type,fieldnum));
        elseif strcmp(set_name,'lauer')
            load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/trilegal/%s_data/%s/isltrilegal_%02d.mat',set_name,tri_type,fieldnum));
        elseif strcmp(set_name,'ghost')
            disp('ERROR: ',set_name,' has not been implemented. Crashing')
            quit cancel %quit code fast
        end
        
        ntri = numel(V);
        
        isltot = zeros(ntri,1);
        islmasked = zeros(ntri,1);
        trimag = cell([ntri,1]); %preallocate cell
        triCoords = cell([ntri,1]); %preallocate cell
        lIltot = cell([ntri,1]); %preallocate cell
        %     lIlcat = cell([ntri,1]); %preallocate cell
        
        for jfile=1:ntri
            
            if strcmp(tri_type,'ubvri')
                
                trimag{jfile} = V(jfile).mlf;
            
                % step 3: convert from LORRI-band mag tmaskso flux
                Fcat = data.cal.vzero .* 10.^(-V(jfile).mlf/2.5);
                
            elseif strcmp(tri_type,'gaia')
                
                trimag{jfile} = V(jfile).G;
            
                % step 3: convert from LORRI-band mag to flux
                Fcat = data.cal.vzero .* 10.^(-V(jfile).G/2.5);
                
            end
            
            %         % step 4: make a mask function
            %         if tri_gaia == 0
            %             whpl = V(jfile).V > max_mag;
            %         elseif tri_gaia == 1
            %             whpl = V(jfile).V > tri_mag;
            %         end
            %
            %         magcat = mag(whpl);
            
            % step 5: convert to surface brightness
            lIltot{jfile} = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea.*(pi./180).^2);%.*(data.astrom.imagew.*data.astrom.imageh);
            %         lIltot{jfile} = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea .* (pi./180).^2);
            %         lIlcat = 1e-26.*1e9.*data.cal.nu.*Fcat(whpl)./(surveyarea .* (pi./180).^2);
            
            isltot(jfile) = sum(lIltot{jfile});
            %         islmasked(jfile) = sum(lIlcat);
            
            % step 6: generate random coordinates for sources
            triCoords{jfile} = rand([length(trimag{jfile}),2])*(inputMap_size(1)-1)+1; %!!!place some sources randomly, they will overlap/be too close to edge [goes from 1 to 256 - for low-res grid]
            
        end
        
        %     trimag = V(1).mlf;
        %
        %     % step 3: convert from LORRI-band mag to flux
        %     Fcat = data.cal.vzero .* 10.^(-V(1).mlf/2.5);
        %
        %     % step 4: convert to surface brightness
        %     lIltot = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea.*(pi./180).^2).*(data.astrom.imagew.*data.astrom.imageh);
        %
        %     isltot = sum(lIltot);
        %
        %     % step 5: generate random coordinates for sources
        %     triCoords = rand([length(trimag),2])*(inputMap_size(1)-1)+1; %!!!place some sources randomly, they will overlap/be too close to edge [goes from 1 to 256 - for low-res grid]
        
        
        %% Make simulated images using catalog or Trilegal sources
        magBinsCat = horzcat((min(floor(mymag)):magBinsStepSize:max(ceil(mymag))-magBinsStepSize)',(min(floor(mymag))+magBinsStepSize:magBinsStepSize:max(ceil(mymag)))'); %get the mag bins spaced 1 apart
        magBinsSizeCat = size(magBinsCat); %get size of mag bins since can't compound calls in matlab
        magBinsTri = cell([ntri,1]); %preallocate cell
        magBinsSizeTri = cell([ntri,1]); %preallocate cell
        for ilegal = 1:ntri
            magBinsTri{ilegal} = horzcat((min(floor(trimag{ilegal})):magBinsStepSize:max(ceil(trimag{ilegal}))-magBinsStepSize)',(min(floor(trimag{ilegal}))+magBinsStepSize:magBinsStepSize:max(ceil(trimag{ilegal})))'); %get the mag bins spaced 1 apart
            magBinsSizeTri{ilegal} = size(magBinsTri{ilegal}); %get size of mag bins since can't compound calls in matlab
        end
        inputMapHolder{ifield}.cat_image = zeros([inputMap_size,magBinsSizeCat(1)]); %preallocate inputMap holder for catalog image
        inputMapHolder{ifield}.bins_cat = magBinsCat; %record the associated mag bins
        inputMapHolder{ifield}.bins_tri = magBinsTri; %record the associated mag bins
        inputMapHolder{ifield}.catmag = mymag; %record the associated star mags
        inputMapHolder{ifield}.trimag = trimag; %record the associated star mags
        inputMapHolder{ifield}.catflux = myflux; %record the associated star fluxes
        inputMapHolder{ifield}.triflux = lIltot; %record the associated star fluxes
        inputMapHolder{ifield}.catypix = myypix; %record the associated star x pix
        inputMapHolder{ifield}.catxpix = myxpix; %record the associated star y pix
        inputMapHolder{ifield}.tripix = triCoords; %record the associated star x and y pix (because there may be several trilegal models, so splitting should be done as needed)
        inputMapHolder{ifield}.catradius = myradius; %record the associated radius
        inputMapHolder{ifield}.filename = datafiles(ifile).name; %record the associated file name
        inputMapHolder{ifield}.fieldnum = goodfields(ifield); %record the associated field number
        %     inputMapHolder{ifield}.starmask = starmask; %record starmask just in case (made with radius + an alg)
        
        % Make catalog image
        for irate = 1:magBinsSizeCat(1)
            
            k = (mymag >= magBinsCat(irate,1)) & (mymag < magBinsCat(irate,2)); %get where magnitudes are within the bin
            curr_mymag = mymag(k); %get only values related to the mag range
            curr_myflux = myflux(k); %get only values related to the mag range
            curr_myypix = myypix(k); %get only values related to the mag range
            curr_myxpix = myxpix(k); %get only values related to the mag range
            %     starNum = 100; %!!!get real num of stars, put here
            %     starFlux = rand([starNum,1])*1000; %!!!random star fluxes here
            %     starCoords = rand([starNum,2])*(inputMap_size(1)-1)+1; %!!!place some sources randomly, they will overlap/be too close to edge [goes from 1 to 256 - for low-res grid]
            %         starNum = length(mymag); %num of stars from data
            %         starFlux = data.cal.sbconv .* myflux; %star fluxes from data
            %         starCoords = horzcat(myypix,myxpix); %low-res coordinates for stars from data
            starNum = length(curr_mymag); %num of stars from data
            inputMapHolder{ifield}.catnum(irate) = starNum; %record the number of stars in bin
            inputMapHolder{ifield}.catsum(irate) = sum(curr_myflux); %record the total flux of stars in bin
            starFlux = curr_myflux; %star fluxes from data
            starCoords = horzcat(curr_myypix,curr_myxpix); %low-res coordinates for stars from data
            
            %just coord calcs below
            starCoords_up = starCoords*factor; %get coords in upscaled coord space
            starCoords_upRound = round(starCoords_up); %round the coords, used for placing sources
            starCoords_upDelta = starCoords_up - starCoords_upRound; %get the deltas, used for sub-pixel shift with RegridderShift
            
            %---make inputMap---
            inputMap_up = zeros(inputMap_upSize); %preallocate inputMap
            
            tempInputCoords = zeros([2,2]); %temp variable that deals with edge cases when placing sources in inputMap [for inputMap]
            tempPsfCoords = zeros([2,2]); %temp variable that deals with edge cases when placing the psf into inputMap [for psf]
            tempPsfSize = zeros([2,2]); %temp variable that has the psf's 1/2 sizes post-shift [so a 50x50 psf shifted 0.5 to right will be a 50x51 psf]
            for i = 1:starNum
                %--shift psf by subpixel amount in starCoords_upDelta--
                %don't have that capability yet
                tempPSF = RegridderShift(psf,starCoords_upDelta(i,:),0,0); %this will be the RegridderShift function called on psf
                sizeTempPSF = size(tempPSF); %size of temp psf that was regridded
                %--get current psf size--
                %x axis delta
                if( starCoords_upDelta(i,1) > 0 )
                    tempPsfSize(1,1) = psfSizeD2(1,1); %regular size
                    tempPsfSize(1,2) = psfSizeD2(1,2)+1; %got bigger here
                elseif( starCoords_upDelta(i,1) < 0 )
                    tempPsfSize(1,1) = psfSizeD2(1,1)+1; %got bigger here
                    tempPsfSize(1,2) = psfSizeD2(1,2); %regular size
                else
                    tempPsfSize(1,1) = psfSizeD2(1,1); %regular size
                    tempPsfSize(1,2) = psfSizeD2(1,2); %regular size
                end
                %y axis delta
                if( starCoords_upDelta(i,2) > 0 )
                    tempPsfSize(2,1) = psfSizeD2(2,1); %regular size
                    tempPsfSize(2,2) = psfSizeD2(2,2)+1; %got bigger here
                elseif( starCoords_upDelta(i,2) < 0 )
                    tempPsfSize(2,1) = psfSizeD2(2,1)+1; %got bigger here
                    tempPsfSize(2,2) = psfSizeD2(2,2); %regular size
                else
                    tempPsfSize(2,1) = psfSizeD2(2,1); %regular size
                    tempPsfSize(2,2) = psfSizeD2(2,2); %regular size
                end
                
                %--Deal with stars on the edges--
                %Check 0 edge, x axis
                if( (starCoords_upRound(i,1) - tempPsfSize(1,1)) <= 0 )
                    tempInputCoords(1,1) = starCoords_upRound(i,1)-1; %make it so it's on the edge at index 1
                    tempPsfCoords(1,1) = 1+tempPsfSize(1,1)-(starCoords_upRound(i,1)-1); %gotta truncate psf to match
                else
                    tempInputCoords(1,1) = tempPsfSize(1,1); %can place psf normally
                    tempPsfCoords(1,1) = 1; %can place psf normally
                end
                %Check -1 edge, x axis
                if( (starCoords_upRound(i,1) + tempPsfSize(1,2)) > inputMap_upSize(1) )
                    tempInputCoords(1,2) = inputMap_upSize(1) - starCoords_upRound(i,1); %[-1 b/c matlab implicitly includes final one]
                    tempPsfCoords(1,2) = sizeTempPSF(1)-(tempPsfSize(1,2)-1-tempInputCoords(1,2)); %gotta truncate psf to match
                else
                    tempInputCoords(1,2) = tempPsfSize(1,2)-1; %can place psf normally [-1 b/c matlab implicitly includes final one]
                    tempPsfCoords(1,2) = sizeTempPSF(1); %can place psf normally
                end
                %Check 0 edge, y axis
                if( (starCoords_upRound(i,2) - tempPsfSize(2,1)) <= 0 )
                    tempInputCoords(2,1) = starCoords_upRound(i,2)-1; %make it so it's on the edge at index 1
                    tempPsfCoords(2,1) = 1+tempPsfSize(2,1)-(starCoords_upRound(i,2)-1); %gotta truncate psf to match
                else
                    tempInputCoords(2,1) = tempPsfSize(2,1); %can place psf normally
                    tempPsfCoords(2,1) = 1; %can place psf normally
                end
                %Check -1 edge, y axis
                if( (starCoords_upRound(i,2) + tempPsfSize(2,2)) > inputMap_upSize(2) )
                    tempInputCoords(2,2) = inputMap_upSize(2) - starCoords_upRound(i,2); %[-1 b/c matlab implicitly includes final one]
                    tempPsfCoords(2,2) = sizeTempPSF(2)-(tempPsfSize(2,2)-1-tempInputCoords(2,2)); %gotta truncate psf to match
                else
                    tempInputCoords(2,2) = tempPsfSize(2,2)-1; %can place psf normally [-1 b/c matlab implicitly includes final one]
                    tempPsfCoords(2,2) = sizeTempPSF(2); %can place psf normally
                end
                
                %finally place the star's psf in input map
                inputMap_up(starCoords_upRound(i,1)-tempInputCoords(1,1):starCoords_upRound(i,1)+tempInputCoords(1,2), ...
                    starCoords_upRound(i,2)-tempInputCoords(2,1):starCoords_upRound(i,2)+tempInputCoords(2,2)) = ...
                    inputMap_up(starCoords_upRound(i,1)-tempInputCoords(1,1):starCoords_upRound(i,1)+tempInputCoords(1,2), ...
                    starCoords_upRound(i,2)-tempInputCoords(2,1):starCoords_upRound(i,2)+tempInputCoords(2,2)) + ...
                    tempPSF(tempPsfCoords(1,1):tempPsfCoords(1,2),tempPsfCoords(2,1):tempPsfCoords(2,2))*starFlux(i); %place that psf, fits it as possible based on edges
            end
            
            inputMap = RegridderZen(inputMap_up,inputMap_size); %downsize the inputMap to the desired size
            
            inputMapHolder{ifield}.cat_image(:,:,irate) = inputMap; %record inputMap
        end
        
        %preallocate for trilegal image
        ilegal_toUse = 1; %make length(trimag) to do all of the trilegal instances
        inputMapHolder{ifield}.tri_image = cell([ilegal_toUse,1]); %preallocate cell
        for ilegal = 1:ilegal_toUse
            
            inputMapHolder{ifield}.tri_image{ilegal} = zeros([inputMap_size,magBinsSizeTri{ilegal}(1)]); %preallocate inputMap holder for trilegal image
        end
        inputMapHolder{ifield}.trinum = cell([length(trimag) ,1]); %preallocate cell [these are needed even if not making all trilegal images]
        inputMapHolder{ifield}.trisum = cell([length(trimag) ,1]); %preallocate cell
        for ilegal = 1:length(trimag)
            inputMapHolder{ifield}.trinum{ilegal} = zeros(magBinsSizeTri{ilegal}(1),1); %preallocate
            inputMapHolder{ifield}.trisum{ilegal} = zeros(magBinsSizeTri{ilegal}(1),1); %preallocate
        end
        
        %Calc important values needed even if not making all trilegal images
        for ilegal = 1:length(trimag)  %change above to do all of the trilegal instances
            for irate = 1:magBinsSizeTri{ilegal}(1)
                k = (trimag{ilegal} >= magBinsTri{ilegal}(irate,1)) & (trimag{ilegal} < magBinsTri{ilegal}(irate,2)); %get where magnitudes are within the bin
                curr_mymag = trimag{ilegal}(k); %get only values related to the mag range
                curr_myflux = lIltot{ilegal}(k); %get only values related to the mag range
                starNum = length(curr_mymag); %num of stars from data
                inputMapHolder{ifield}.trinum{ilegal}(irate) = starNum; %record the number of stars in bin
                if( sum(k) == 0 )
                    inputMapHolder{ifield}.trisum{ilegal}(irate) = NaN;
                else
                    inputMapHolder{ifield}.trisum{ilegal}(irate) = sum(curr_myflux); %record the total flux of stars in bin
                end
            end
        end
        
        %Make trilegal image
        for ilegal = 1:ilegal_toUse %change above to do all of the trilegal instances
            for irate = 1:magBinsSizeTri{ilegal}(1)
                
                k = (trimag{ilegal} >= magBinsTri{ilegal}(irate,1)) & (trimag{ilegal} < magBinsTri{ilegal}(irate,2)); %get where magnitudes are within the bin
                curr_mymag = trimag{ilegal}(k); %get only values related to the mag range
                curr_myflux = lIltot{ilegal}(k); %get only values related to the mag range
                curr_myypix = triCoords{ilegal}(k,1); %get only values related to the mag range
                curr_myxpix = triCoords{ilegal}(k,2); %get only values related to the mag range
                %     starNum = 100; %!!!get real num of stars, put here
                %     starFlux = rand([starNum,1])*1000; %!!!random star fluxes here
                %     starCoords = rand([starNum,2])*(inputMap_size(1)-1)+1; %!!!place some sources randomly, they will overlap/be too close to edge [goes from 1 to 256 - for low-res grid]
                %         starNum = length(mymag); %num of stars from data
                %         starFlux = data.cal.sbconv .* myflux; %star fluxes from data
                %         starCoords = horzcat(myypix,myxpix); %low-res coordinates for stars from data
                starNum = length(curr_mymag); %num of stars from data
                %             inputMapHolder{ifield}.trinum{ilegal}(irate) = starNum; %record the number of stars in bin
                %             inputMapHolder{ifield}.trisum{ilegal}(irate) = sum(curr_myflux); %record the total flux of stars in bin
                starFlux = curr_myflux; %star fluxes from data
                starCoords = horzcat(curr_myypix,curr_myxpix); %low-res coordinates for stars from data
                
                %just coord calcs below
                starCoords_up = starCoords*factor; %get coords in upscaled coord space
                starCoords_upRound = round(starCoords_up); %round the coords, used for placing sources
                starCoords_upDelta = starCoords_up - starCoords_upRound; %get the deltas, used for sub-pixel shift with RegridderShift
                
                %---make inputMap---
                inputMap_up = zeros(inputMap_upSize); %preallocate inputMap
                
                tempInputCoords = zeros([2,2]); %temp variable that deals with edge cases when placing sources in inputMap [for inputMap]
                tempPsfCoords = zeros([2,2]); %temp variable that deals with edge cases when placing the psf into inputMap [for psf]
                tempPsfSize = zeros([2,2]); %temp variable that has the psf's 1/2 sizes post-shift [so a 50x50 psf shifted 0.5 to right will be a 50x51 psf]
                for i = 1:starNum
                    %--shift psf by subpixel amount in starCoords_upDelta--
                    %don't have that capability yet
                    tempPSF = RegridderShift(psf,starCoords_upDelta(i,:),0,0); %this will be the RegridderShift function called on psf
                    sizeTempPSF = size(tempPSF); %size of temp psf that was regridded
                    %--get current psf size--
                    %x axis delta
                    if( starCoords_upDelta(i,1) > 0 )
                        tempPsfSize(1,1) = psfSizeD2(1,1); %regular size
                        tempPsfSize(1,2) = psfSizeD2(1,2)+1; %got bigger here
                    elseif( starCoords_upDelta(i,1) < 0 )
                        tempPsfSize(1,1) = psfSizeD2(1,1)+1; %got bigger here
                        tempPsfSize(1,2) = psfSizeD2(1,2); %regular size
                    else
                        tempPsfSize(1,1) = psfSizeD2(1,1); %regular size
                        tempPsfSize(1,2) = psfSizeD2(1,2); %regular size
                    end
                    %y axis delta
                    if( starCoords_upDelta(i,2) > 0 )
                        tempPsfSize(2,1) = psfSizeD2(2,1); %regular size
                        tempPsfSize(2,2) = psfSizeD2(2,2)+1; %got bigger here
                    elseif( starCoords_upDelta(i,2) < 0 )
                        tempPsfSize(2,1) = psfSizeD2(2,1)+1; %got bigger here
                        tempPsfSize(2,2) = psfSizeD2(2,2); %regular size
                    else
                        tempPsfSize(2,1) = psfSizeD2(2,1); %regular size
                        tempPsfSize(2,2) = psfSizeD2(2,2); %regular size
                    end
                    
                    %--Deal with stars on the edges--
                    %Check 0 edge, x axis
                    if( (starCoords_upRound(i,1) - tempPsfSize(1,1)) <= 0 )
                        tempInputCoords(1,1) = starCoords_upRound(i,1)-1; %make it so it's on the edge at index 1
                        tempPsfCoords(1,1) = 1+tempPsfSize(1,1)-(starCoords_upRound(i,1)-1); %gotta truncate psf to match
                    else
                        tempInputCoords(1,1) = tempPsfSize(1,1); %can place psf normally
                        tempPsfCoords(1,1) = 1; %can place psf normally
                    end
                    %Check -1 edge, x axis
                    if( (starCoords_upRound(i,1) + tempPsfSize(1,2)) > inputMap_upSize(1) )
                        tempInputCoords(1,2) = inputMap_upSize(1) - starCoords_upRound(i,1); %[-1 b/c matlab implicitly includes final one]
                        tempPsfCoords(1,2) = sizeTempPSF(1)-(tempPsfSize(1,2)-1-tempInputCoords(1,2)); %gotta truncate psf to match
                    else
                        tempInputCoords(1,2) = tempPsfSize(1,2)-1; %can place psf normally [-1 b/c matlab implicitly includes final one]
                        tempPsfCoords(1,2) = sizeTempPSF(1); %can place psf normally
                    end
                    %Check 0 edge, y axis
                    if( (starCoords_upRound(i,2) - tempPsfSize(2,1)) <= 0 )
                        tempInputCoords(2,1) = starCoords_upRound(i,2)-1; %make it so it's on the edge at index 1
                        tempPsfCoords(2,1) = 1+tempPsfSize(2,1)-(starCoords_upRound(i,2)-1); %gotta truncate psf to match
                    else
                        tempInputCoords(2,1) = tempPsfSize(2,1); %can place psf normally
                        tempPsfCoords(2,1) = 1; %can place psf normally
                    end
                    %Check -1 edge, y axis
                    if( (starCoords_upRound(i,2) + tempPsfSize(2,2)) > inputMap_upSize(2) )
                        tempInputCoords(2,2) = inputMap_upSize(2) - starCoords_upRound(i,2); %[-1 b/c matlab implicitly includes final one]
                        tempPsfCoords(2,2) = sizeTempPSF(2)-(tempPsfSize(2,2)-1-tempInputCoords(2,2)); %gotta truncate psf to match
                    else
                        tempInputCoords(2,2) = tempPsfSize(2,2)-1; %can place psf normally [-1 b/c matlab implicitly includes final one]
                        tempPsfCoords(2,2) = sizeTempPSF(2); %can place psf normally
                    end
                    
                    %finally place the star's psf in input map
                    inputMap_up(starCoords_upRound(i,1)-tempInputCoords(1,1):starCoords_upRound(i,1)+tempInputCoords(1,2), ...
                        starCoords_upRound(i,2)-tempInputCoords(2,1):starCoords_upRound(i,2)+tempInputCoords(2,2)) = ...
                        inputMap_up(starCoords_upRound(i,1)-tempInputCoords(1,1):starCoords_upRound(i,1)+tempInputCoords(1,2), ...
                        starCoords_upRound(i,2)-tempInputCoords(2,1):starCoords_upRound(i,2)+tempInputCoords(2,2)) + ...
                        tempPSF(tempPsfCoords(1,1):tempPsfCoords(1,2),tempPsfCoords(2,1):tempPsfCoords(2,2))*starFlux(i); %place that psf, fits it as possible based on edges
                end
                
                inputMap = RegridderZen(inputMap_up,inputMap_size); %downsize the inputMap to the desired size
                
                inputMapHolder{ifield}.tri_image{ilegal}(:,:,irate) = inputMap; %record inputMap
            end
        end
        
        %% Make ISL plots

        % Plot integrated surface brightness from various sources
        if FLG_firstRun
            h = figure();
            clf;
            hold on
        end
        % legend_handle(end+1) = scatter(inputMapHolder{ifield}.bins_cat(:,1)+0.5,inputMapHolder{ifield}.catsum,70,'filled',...
        %     'MarkerEdgeColor',PLOT_color{colorCntr},'MarkerFaceColor',PLOT_color{colorCntr}); %plot summed flux
        if( FLG_cumsum == 0 )
            legend_handle(end+1) = plot(inputMapHolder{ifield}.bins_cat(:,1)+0.5,inputMapHolder{ifield}.catsum,...
                'color', PLOT_color{colorCntr}, 'LineWidth', 2); %plot summed flux
        else
            legend_handle(end+1) = plot(inputMapHolder{ifield}.bins_cat(:,1)+0.5,flip(cumsum(flip(inputMapHolder{ifield}.catsum))),...
                'color', PLOT_color{colorCntr}, 'LineWidth', 2); %plot summed flux
        end
%         scatter(inputMapHolder{ifield}.bins_cat(:,1)+0.5,inputMapHolder{ifield}.catnum,70,'filled'); %plot number counts
        binsMax = inputMapHolder{ifield}.bins_tri{1}(:,1); %start off binsMax
        for ilegal = 1:length(trimag) %do all of the trilegal instances
            binsTemp = inputMapHolder{ifield}.bins_tri{ilegal}(:,1); %get this as temp so can slice n dice it
            binsLogic = ismember(binsTemp,binsMax); %get if 1st array is contained by 2nd
            if( all(binsLogic) == 0 )
                binsMax = vertcat(binsMax,binsTemp(~binsLogic)); %get the maximum bins that can be found in all of the trilegal models
            end
        end
        binsMax = sort(binsMax); %sort because might not be in right order
        trisum_comb = zeros(length(trimag),length(binsMax));
        trinum_comb = zeros(length(trimag),length(binsMax));
        for ilegal = 1:length(trimag) %do all of the trilegal instances
            bins_temp = inputMapHolder{ifield}.bins_tri{ilegal}(:,1);
            trisum_temp = inputMapHolder{ifield}.trisum{ilegal};
            trinum_temp = inputMapHolder{ifield}.trinum{ilegal};
            index_start = find( bins_temp(1) == binsMax); %the trilegal bins always end at 33-34 bin, but can start at different magnitudes
            trisum_temp = padarray(trisum_temp,index_start-1,NaN,'pre'); %record, pad NaNs
            trinum_temp = padarray(trinum_temp,index_start-1,NaN,'pre'); %record, pad NaNs
            index_end = find( bins_temp(end) == binsMax ); %the trilegal bins always end at 33-34 bin, but can start at different magnitudes
            trisum_comb(ilegal,:) = padarray(trisum_temp,length(binsMax)-index_end,NaN,'post'); %record, pad NaNs
            trinum_comb(ilegal,:) = padarray(trinum_temp,length(binsMax)-index_end,NaN,'post'); %record, pad NaNs
        end
%             errorbar(binsMax+0.5,nanmean(trisum_comb,1)',nanstd(trisum_comb,1)','o'); %mean/std of summed flux
        % legend_handle(end+1) = errorbar(binsMax+0.5,nanmean(trisum_comb,1)',nanmin(trisum_comb)',nanmax(trisum_comb)','o',...
        %     'Color',PLOT_color{colorCntr},'MarkerEdgeColor',PLOT_color{colorCntr},'MarkerFaceColor',PLOT_color{colorCntr}); %min/max of summed flux
        
        binsMax_plot = binsMax; %regular
        if( FLG_cumsum == 1 )
            trisum_comb_nanloc = isnan(trisum_comb); %find the nan values to ensure they stay nan
            trisum_comb = flip(cumsum(flip(trisum_comb,2),2,'omitnan'),2); %ignore nans in cumsum
            trisum_comb(trisum_comb_nanloc) = nan; %return any nan'd values to nan
        end
        trisum_comb_mean = nanmean(trisum_comb,1)';
        trisum_comb_max = trisum_comb_mean+nanmax(trisum_comb)';
        trisum_comb_min = trisum_comb_mean-nanmin(trisum_comb)';
        if( any(isnan(nanmean(trisum_comb,1))) )
            jk = ~isnan(nanmean(trisum_comb,1));
            binsMax_plot = binsMax_plot(jk); %regular
            trisum_comb_max = trisum_comb_max(jk);
            trisum_comb_min = trisum_comb_min(jk);
            trisum_comb_mean = trisum_comb_mean(jk); %for later
        end
        trisum_comb_min(trisum_comb_min == 0) = trisum_comb_mean(trisum_comb_min == 0); %if 0 make it the mean (happens if only 1 pt)
        legend_handle(end+1) = fill([binsMax_plot+0.5; flip(binsMax_plot+0.5)], ...
            [trisum_comb_min; flip(trisum_comb_max)],...
            sscanf(PLOT_color{colorCntr}(2:end),'%2x%2x%2x',[1 3])/255, ...
            'FaceColor', sscanf(PLOT_color{colorCntr}(2:end),'%2x%2x%2x',[1 3])/255, 'FaceAlpha', 0.22, 'EdgeColor', 'none'); %, 'EdgeColor', 'none'
%         errorbar(binsMax+0.5,nanmean(trinum_comb,1)',nanstd(trinum_comb)','o'); %mean/std of number counts
        %     scatter(wing_list(:,1),wing_list(:,2),'filled');
        %     scatter(tri_only_list(:,1),tri_only_list(:,2),70,'filled');
        %     scatter(tri_list(:,1),tri_list(:,2),50,'filled');33
        %     scatter(gaia_list(:,1),gaia_list(:,2),30,'filled');
        
        %--legend accumulator---
        if use_gaia == 1
            legend_text{end+1} = ['{\it Gaia} Field ',num2str(uniqueField_id(colorCntr))];
            legend_text{end+1} = ['TRILEGAL Field ',num2str(uniqueField_id(colorCntr))];
            % legend('Gaia Catalog',sprintf('Trilegal %s',tri_type));
        else
            legend_text{end+1} = ['USNOB1 ',set_name,num2str(goodfields(ifield))];
            legend_text{end+1} = ['TRILEGAL ',tri_type,' ',set_name,num2str(goodfields(ifield))];
            % legend('USNOB1 Catalog',sprintf('Trilegal %s',tri_type));
        end

        if FLG_firstRun
            set(gca,'YScale','log'); %use log on y axis
            %     title(sprintf('m = %1.1f',mag));
            xlabel('Mag Bin');
            if( FLG_cumsum == 0 )
                ylabel('Surface Brightness in Bin [nW m^-2 sr^-1]');
            else
                ylabel('Cumulative Surface Brightness in Bin [nW m^-2 sr^-1]');
            end
    %         ylabel('Num. of Stars');
            FLG_firstRun = false; %turn off first run
        end
        colorCntr = colorCntr + 1;
    end %end ifield loop (3/5/6/...)
end %end iset loop (old/new/...)
%reorder so lines ontop of patches
zorder_curr = flip(get(gca, 'Children'));
zorder = gobjects(size(zorder_curr));
zorder_frontCntr = 1; %put it at the front
zorder_backCntr = length(zorder_curr)/2+1; %put it in the middle
for i=1:length(zorder_curr)
    if( strcmp(zorder_curr(i).Type,'patch') )
        zorder(zorder_backCntr) = zorder_curr(i); %move it to the back
        zorder_backCntr = zorder_backCntr + 1; %increment
    else
        zorder(zorder_frontCntr) = zorder_curr(i); %to the front
        zorder_frontCntr = zorder_frontCntr + 1; %increment
    end
end
set(gca, 'Children',zorder)

%finish plotting by applying legend
legend(legend_handle,legend_text)
%save plot maybe then
ext = '.png';
if use_gaia == 1
    imagename = [pwd,'/imagesim_field_',char(join(dataSets.sets,'_')),'_Gaia',ext];
else
    imagename = [pwd,'/imagesim_field_',char(join(dataSets.sets,'_')),'_USNOB1',ext];
end
print(h,imagename, '-dpng');
display('donezo');
