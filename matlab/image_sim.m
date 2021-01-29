%make simulated lorri images for the purpose of estimating flux
%contributions from masking and integrated starlight

clc;
clear variables;
close all;

%set pipeline settings
use_gaia = 1;
max_mag = 30;

magBinsStepSize = 1; %sets the step size for the mag bins

paths = get_paths_new();

datafiles = dir(sprintf('%s*.mat',paths.datadir));

nfiles = size(datafiles,1);

inputMapHolder = cell([nfiles,1]); %preallocate inputMap holder - holds a struct for each file

%% from catalog, retrieve lists of stars with mags and coords
for ifile=1:nfiles
    fprintf(['On file ',num2str(ifile),' of ',num2str(size(datafiles,1)),'.\n']);
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    if data.header.fieldnum > 1
        
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
                    rangeexclude = rangeexclude+1;
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
                    numinimg = numinimg + 1;
                    
                    % require that the magnitude is within sensible bounds
                    if thismag < max_mag & ~isnan(thismag) & thissg > 0
                        numinbnds = numinbnds + 1;
                        
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
        
        % step 2: read in the trilegal information
        load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/isltrilegal_%02d.mat',fieldnum));
        
        ntri = numel(V);
        
        isltot = zeros(ntri,1);
        islmasked = zeros(ntri,1);
        trimag = cell([ntri,1]); %preallocate cell
        triCoords = cell([ntri,1]); %preallocate cell
        lIltot = cell([ntri,1]); %preallocate cell
        %     lIlcat = cell([ntri,1]); %preallocate cell
        
        for jfile=1:ntri
            
            trimag{jfile} = V(jfile).mlf;
            
            % step 3: convert from LORRI-band mag to flux
            Fcat = data.cal.vzero .* 10.^(-V(jfile).mlf/2.5);
            
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
        inputMapHolder{ifile}.cat_image = zeros([inputMap_size,magBinsSizeCat(1)]); %preallocate inputMap holder for catalog image
        inputMapHolder{ifile}.bins_cat = magBinsCat; %record the associated mag bins
        inputMapHolder{ifile}.bins_tri = magBinsTri; %record the associated mag bins
        inputMapHolder{ifile}.catmag = mymag; %record the associated star mags
        inputMapHolder{ifile}.trimag = trimag; %record the associated star mags
        inputMapHolder{ifile}.catflux = myflux; %record the associated star fluxes
        inputMapHolder{ifile}.triflux = lIltot; %record the associated star fluxes
        inputMapHolder{ifile}.catypix = myypix; %record the associated star x pix
        inputMapHolder{ifile}.catxpix = myxpix; %record the associated star y pix
        inputMapHolder{ifile}.tripix = triCoords; %record the associated star x and y pix (because there may be several trilegal models, so splitting should be done as needed)
        inputMapHolder{ifile}.catradius = myradius; %record the associated radius
        inputMapHolder{ifile}.filename = datafiles(ifile).name; %record the associated file name
        %     inputMapHolder{ifile}.starmask = starmask; %record starmask just in case (made with radius + an alg)
        
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
            inputMapHolder{ifile}.catnum(irate) = starNum; %record the number of stars in bin
            inputMapHolder{ifile}.catsum(irate) = sum(curr_myflux); %record the total flux of stars in bin
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
            
            inputMapHolder{ifile}.cat_image(:,:,irate) = inputMap; %record inputMap
        end
        
        %preallocate for trilegal image
        ilegal_toUse = 1; %make length(trimag) to do all of the trilegal instances
        inputMapHolder{ifile}.tri_image = cell([ilegal_toUse,1]); %preallocate cell
        for ilegal = 1:ilegal_toUse
            
            inputMapHolder{ifile}.tri_image{ilegal} = zeros([inputMap_size,magBinsSizeTri{ilegal}(1)]); %preallocate inputMap holder for trilegal image
        end
        inputMapHolder{ifile}.trinum = cell([length(trimag) ,1]); %preallocate cell [these are needed even if not making all trilegal images]
        inputMapHolder{ifile}.trisum = cell([length(trimag) ,1]); %preallocate cell
        for ilegal = 1:length(trimag)
            inputMapHolder{ifile}.trinum{ilegal} = zeros(magBinsSizeTri{ilegal}(1),1); %preallocate
            inputMapHolder{ifile}.trisum{ilegal} = zeros(magBinsSizeTri{ilegal}(1),1); %preallocate
        end
        
        %Calc important values needed even if not making all trilegal images
        for ilegal = 1:length(trimag)  %change above to do all of the trilegal instances
            for irate = 1:magBinsSizeTri{ilegal}(1)
                k = (trimag{ilegal} >= magBinsTri{ilegal}(irate,1)) & (trimag{ilegal} < magBinsTri{ilegal}(irate,2)); %get where magnitudes are within the bin
                curr_mymag = trimag{ilegal}(k); %get only values related to the mag range
                curr_myflux = lIltot{ilegal}(k); %get only values related to the mag range
                starNum = length(curr_mymag); %num of stars from data
                inputMapHolder{ifile}.trinum{ilegal}(irate) = starNum; %record the number of stars in bin
                if( sum(k) == 0 )
                    inputMapHolder{ifile}.trisum{ilegal}(irate) = NaN;
                else
                    inputMapHolder{ifile}.trisum{ilegal}(irate) = sum(curr_myflux); %record the total flux of stars in bin
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
                %             inputMapHolder{ifile}.trinum{ilegal}(irate) = starNum; %record the number of stars in bin
                %             inputMapHolder{ifile}.trisum{ilegal}(irate) = sum(curr_myflux); %record the total flux of stars in bin
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
                
                inputMapHolder{ifile}.tri_image{ilegal}(:,:,irate) = inputMap; %record inputMap
            end
        end
        
        %% Make ISL plots
        
        % Plot integrated surface brightness from various sources
        h = figure(1);
        clf;
        hold on
        scatter(inputMapHolder{ifile}.bins_cat(:,1)+0.5,inputMapHolder{ifile}.catsum,70,'filled'); %plot summed flux
%         scatter(inputMapHolder{ifile}.bins_cat(:,1)+0.5,inputMapHolder{ifile}.catnum,70,'filled'); %plot number counts
        binsMax = inputMapHolder{ifile}.bins_tri{1}(:,1); %start off binsMax
        for ilegal = 1:length(trimag) %do all of the trilegal instances
            binsTemp = inputMapHolder{ifile}.bins_tri{ilegal}(:,1); %get this as temp so can slice n dice it
            binsLogic = ismember(binsTemp,binsMax); %get if 1st array is contained by 2nd
            if( all(binsLogic) == 0 )
                binsMax = vertcat(binsMax,binsTemp(~binsLogic)); %get the maximum bins that can be found in all of the trilegal models
            end
        end
        binsMax = sort(binsMax); %sort because might not be in right order
        trisum_comb = zeros(length(trimag),length(binsMax));
        trinum_comb = zeros(length(trimag),length(binsMax));
        for ilegal = 1:length(trimag) %do all of the trilegal instances
            bins_temp = inputMapHolder{ifile}.bins_tri{ilegal}(:,1);
            trisum_temp = inputMapHolder{ifile}.trisum{ilegal};
            trinum_temp = inputMapHolder{ifile}.trinum{ilegal};
            index_start = find( bins_temp(1) == binsMax); %the trilegal bins always end at 33-34 bin, but can start at different magnitudes
            trisum_temp = padarray(trisum_temp,index_start-1,NaN,'pre'); %record, pad NaNs
            trinum_temp = padarray(trinum_temp,index_start-1,NaN,'pre'); %record, pad NaNs
            index_end = find( bins_temp(end) == binsMax ); %the trilegal bins always end at 33-34 bin, but can start at different magnitudes
            trisum_comb(ilegal,:) = padarray(trisum_temp,length(binsMax)-index_end,NaN,'post'); %record, pad NaNs
            trinum_comb(ilegal,:) = padarray(trinum_temp,length(binsMax)-index_end,NaN,'post'); %record, pad NaNs
        end
%             errorbar(binsMax+0.5,nanmean(trisum_comb,1)',nanstd(trisum_comb,1)','o'); %mean/std of summed flux
        errorbar(binsMax+0.5,nanmean(trisum_comb,1)',nanmin(trisum_comb)',nanmax(trisum_comb)','o'); %min/max of summed flux
%         errorbar(binsMax+0.5,nanmean(trinum_comb,1)',nanstd(trinum_comb)','o'); %mean/std of number counts
        %     scatter(wing_list(:,1),wing_list(:,2),'filled');
        %     scatter(tri_only_list(:,1),tri_only_list(:,2),70,'filled');
        %     scatter(tri_list(:,1),tri_list(:,2),50,'filled');
        %     scatter(gaia_list(:,1),gaia_list(:,2),30,'filled');
        set(gca,'YScale','log'); %use log on y axis
        legend('Gaia Catalog','Trilegal');
        %     title(sprintf('m = %1.1f',mag));
        xlabel('Mag Bin');
        ylabel('Surface Brightness in Bin [nW m^-2 sr^-1]');
%         ylabel('Num. of Stars');
        ext = '.png';
        imagename = sprintf('%s%s%s',paths.simdir,data.header.timestamp,ext);
        print(h,imagename, '-dpng');
        display('donezo');
    end
end