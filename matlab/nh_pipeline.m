%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_pipeline.m
%%  Jun 2016, MZ
%%  For each file in the NH data directories, this program:
%%   1) Reads it in.
%%   2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nh_pipeline()
clear all
close all
%get paths for new data files or old data files
% paths = get_paths_old_ghosts();
% paths = get_paths_old();
% paths = get_paths_new();
% paths = get_paths_lauer();
paths = get_paths_newest();

datadir = paths.datadir;
% Save which data we're looking at
if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/ghosts/') == 1
    data_type = 'ghost';
    start = 1;
elseif strcmp(paths.datadir,'/data/symons/NH_old_data/mat/good/') == 1
    data_type = 'old';
    start = 1;
elseif strcmp(paths.datadir,'/data/symons/nh_data/mat/') == 1
    data_type = 'new';
    start = 40; % Skip first 40 images (field 1)
elseif strcmp(paths.datadir,'/data/symons/nh_data_lauer/mat/') == 1 || strcmp(paths.datadir,'/data/symons/nh_data_new/mat/') 
    data_type = 'new';
    start = 1;
end

%% Set user parameters here
%choose method (old, old_corr, new)
flag_method = 'new';

%set portion of pipeline to run
procstepflag = 1; %1 - masking, 2 - add meta, 3 - jail bars, 4 - calibrate, 5 - diff ghosts & scattering, 6 - isl, 7 - dgl, 8 - extinction

% Set error flags - error on means error of this type is being added in
errflag_mags = 0; %1 is on, 0 is off - run 1-5
errflag_psf = 0; % - run 6
errflag_gals = 0; % - run 1-4
% Set number of MC iterations for gal error
if errflag_gals == 1
    mc_it = 100;
else
    mc_it = 1;
end

% Set error string, used to create data substruct
if errflag_mags == 1
    err_str = 'err_mags';
    err_on = 1;
elseif errflag_psf == 1
    err_str = 'err_psf';
    err_on = 1;
elseif errflag_gals == 1
    err_str = 'err_gals';
    err_on = 1;
else
    err_str = 'no_err';
    err_on = 0;
end

% Set and save run parameters to params struct
params.data_type = data_type;
params.method = flag_method;
params.step = procstepflag;
params.err_on = err_on;
params.err_mags = errflag_mags;
params.err_psf = errflag_psf;
params.err_gals = errflag_gals;
params.err_str = err_str;

% Save run parameters
save('run_params.mat','params');

%use USNOB1 or Gaia for star masking
if strcmp(flag_method,'new') == 1
    use_gaia = 1; %0 for USNOB1
elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'old') == 1)
    use_gaia = 0;
end

% If method file does not exist (data has not been run before), set
% old_method to 'none' so star mask is created
if( ~(isfile(sprintf('%smethod.txt',data_type))))
    old_method = 'none';
else
    % If method file exists, read saved text file for data_type and see which method last used
    fileID = fopen([data_type,'method.txt'],'r');
    old_method = fscanf(fileID,'%s');
    fclose(fileID);
end

% If method has changed, redo star mask
if strcmp(flag_method,old_method) == 1
    new_star_mask = 0; %0 for skipping star masking and using existing mask
end
if strcmp(flag_method,old_method) == 0 || errflag_mags == 1 || errflag_gals == 1
    new_star_mask = 1;
end
% new_star_mask = 1; %!!!ALWAYS MAKE NEW STAR MASK - TURN OFF FOR BIG RUN!!!

% Save text file per data_type listing most recently run method
fileID = fopen([data_type,'method.txt'],'w');
fprintf(fileID,'%s',flag_method);
fclose(fileID);

%use Gaia or just Trilegal for ISL calc
tri_gaia = 0; %0 for only Trilegal, 1 for Gaia and Trilegal - only need Trilegal if using Gaia for masking

% Choose type of Trilegal magnitudes to use
tri_type = 'gaia'; % Select 'ubvri' for interpolated LORRI mags or 'gaia' for G mags

%keep or omit galaxies from Gaia data
gals = 0; %0 for galaxies removed, 1 for galaxies included - galaxies removed is fine, little effect

%set max mag for masking (mask all stars brighter than this), also
%calculate gaia isl between this and 20
if strcmp(flag_method,'new') == 1
    max_mag = 21; %17.75 old value, 21 new value
elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'old') == 1)
    max_mag = 17.75; % if old method, this is overwritten in makemask
end

%set min mag for trilegal isl (calculate isl for stars fainter than this)
tri_mag = 20; % Only used when tri_gaia == 1 (Gaia and Trilegal used for ISL)

%set max mag for psf wing isl calc (from USNOB1 or Gaia based on masking catalog)
wing_mag = max_mag;

%optionally save lists of stars and flux
save_file = 0; %1 to save


%%
%load directory of saved mat files for all images
datafiles = dir(sprintf('%s*.mat',datadir));

mydate = zeros(size(datafiles));
mytemp = zeros(size(datafiles));
mymean = zeros(size(datafiles));
myref = zeros(size(datafiles));
myeng = zeros(size(datafiles));
myisl = zeros(size(datafiles));

photcurr = zeros(size(datafiles));
checkmax = zeros(size(datafiles));
checkmin = zeros(size(datafiles));

ghostdist = zeros(size(datafiles));

totalghosts = 0;
totalrealghosts = 0;

if procstepflag == 3
    evenMinusOdd = zeros( size(datafiles,1), 1); %preallocate
    evenMinusOddFixed = zeros( size(datafiles,1), 1); %preallocate
    oddMinusEvenref = zeros( size(datafiles,1), 1);
    oddMinusEvenrefFixed = zeros( size(datafiles,1), 1);
end

%loop through all images
timePerSum = 0; %for time estimating
timeCntr = 0;
for ifile=start:size(datafiles,1)
    tic
    fprintf('\nOn file %d of %d\n',ifile,size(datafiles,1));
    if( timeCntr ~= 0 )
        timePerSum = timePerSum + timePer;
        fprintf('Time per file: %f sec, ETA: %f min\n',timePerSum/(timeCntr),timePerSum/(timeCntr)*(size(datafiles,1)-ifile-1)/60);
    end

    %load individual mat file for image
    load(sprintf('%s%s',datadir,datafiles(ifile).name));
    % Display file name
    disp(data.header.rawfile)

    % Display field number
    disp(data.header.fieldnum)

    %ensure error struct is made if error is on
    if( (params.err_on == 1) && (isfield(data,params.err_str) == 0) )
        data.(params.err_str) = struct;
    end

    % Stop on a certain file
    if ifile == 41  %55
        %     if ifile == 65 | ifile == 139 | ifile == 252 | ifile == 281 | ifile == 341
        fprintf('ahhhh')
        %         data.header.bad = 1;
    else
        %         data.header.bad = 0;
    end

    % Stop on a certain field
    %     if data.header.fieldnum == 3
    %         fprintf('ahhhh')
    %     end

    % Analyze just one specific file
%     if strcmp(data.header.timestamp,'2458731.2935809') == 1

    % Restrict newest data to only good fields
    if (strcmp(paths.datadir,'/data/symons/nh_data_new/mat/') == 1) && any(data.header.fieldnum == [2,4,5,6,7,12,15,16,17,19,20,22,23])
    
    % If running err_gals MC, do more than one iteration
    for mc = 1:mc_it
    fprintf('Current MC num: %i of %i | On file %d of %d\n',mc,mc_it,ifile,size(datafiles,1));
        if procstepflag >= 1
            fprintf('Masking and ghost analysis\n');
            %% Masking and optical ghosts

            %remove recursive err_str
            if errflag_mags == 1
                if( isfield(data.(params.err_str),params.err_str) == 1 )
                    data.(params.err_str) = rmfield(data.(params.err_str),params.err_str); %remove recursion
                end
            end

            if strcmp(flag_method,'new') == 1 %!!!NEED TO TURN BACK ON!!!
                %save ra, dec, pix, and mag of stars in image that could be causing ghosts
                [data, ghostcount, realghostcount] = nh_findghoststar(data,paths,params,use_gaia,errflag_mags,mc);
%                 totalghosts = totalghosts + ghostcount;
%                 totalrealghosts = totalrealghosts + realghostcount;
            end

            %
            %manually mask portions of image - not needed for old ghost files
            %         data = nh_make_manmask(data,paths);

            %mask out stars and ghosts in image (also stat and clip mask)
            %may need to redo catalog depending on gals
            %         catalog_data_gaia(gals,paths)
            %         catalog_data_gaia_wide(gals,paths)
            data = nh_makemask(data,paths,params,3,use_gaia,new_star_mask, max_mag, save_file, flag_method, errflag_mags, mc);

            %         maskim = data.data.*~data.mask.onemask;
            %         maskim(maskim==0) = NaN;
            %         checkmax(ifile) = nanmax(nanmax(maskim));
            %         checkmin(ifile) = nanmin(nanmin(maskim));

            %print total masked surface brightness (sum of all masked pixels)
            %         data.cal.sbconv .* sum(sum(data.data(data.mask.onemask)./ data.astrom.exptime))

            %print mask fraction
            %         data.mask.maskfrac

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end

        if procstepflag >= 2
            fprintf('Meta data\n');
            %% Meta data
            %write header values, constants, coordinates, and dark reference pixel data to data file
            %(uses mask, may need to redo mask first)
            %changes values used in calibration
            data = nh_add_meta(data,flag_method); %!!!NEED TO TURN THIS BACK ON !!!

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end

        if procstepflag >= 3
            fprintf('Jail bars\n');
            %% Jail bar correction
            %perform jail bar correction on data.data and data.ref - need mask first, but can redo from original data
            [data, evenMinusOdd(ifile), evenMinusOddFixed(ifile), oddMinusEvenref(ifile), oddMinusEvenrefFixed(ifile)] = nh_jail_bars(data,paths,params,flag_method);

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end

        if procstepflag >= 4
            fprintf('Calibration\n');
            %% Calibration

            %save calibrated image in surface brightness units (nh_calcref saves refcorr which is used here)
            %nh_calcref must be run first *FOR OLD METHOD ONLY*
            data = nh_calibrate(data,paths,params,flag_method,mc);

            if ifile == 1
                % Save conversion factor to text file to be read in by python
                fileID = fopen('conv.txt','w');
                fprintf(fileID,'%.16f',data.cal.sbconv);
                fclose(fileID);
            end

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end
        
        if procstepflag == 5
            fprintf('Diff ghosts and scattering\n');
            %% Diffuse ghosts and scattering

            if strcmp(flag_method,'new') == 1
                %Calculate diffuse contribution from all stars in range to cause a
                %ghost. List of stars from nh_findghoststar. Later subtracted from
                %image mean.
                data = nh_diffghost(data,paths,params); %!!!NEED TO TURN THIS BACK ON !!!
                %Calculate extended diffuse scattering from all-sky ISL
                %Depends on calibration, redo if calibration changes
                if errflag_mags ~= 1
                    data = nh_scattering(data, paths, errflag_mags, params);
                end
            end

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end

        if procstepflag == 6
            fprintf('ISL calculation\n');
            %% ISL calculation

            %calculate ISL from USNOB1 and Trilegal - depends on mask, need to redo if new mask
            %may need to redo catalog depending on gals
            %         catalog_data_gaia(gals,paths)
            % Need to run nh_stack_psf() first to save psf to data
            data = nh_calcisl(data, paths, params, use_gaia, tri_gaia, tri_mag, wing_mag, max_mag, save_file, flag_method, errflag_psf, tri_type);

            %         wing = data.isl.usnowing
            %         triisl = data.isl.trimeanmasksize
            %         gaiaisl = data.isl.gaiamean

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end

        if procstepflag == 7
            fprintf('DGL calculation\n');
            %% DGL calculation
            %         if data.header.fieldnum == 8
            %             fprintf('its here');
            %         end

            %calculate DGL using Planck (or IRIS) data
            data = nh_calcdgl(data, paths, flag_method);

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end

        if procstepflag == 8
            fprintf('Extinction calculation\n');
            %% Extinction calculation

            %calculate extinction using SFD data
            data = nh_extinction(data, paths);

            %overwrite data file with changes
            save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');

        end
    
    end % END FOR MC it

% end % END IF running one file
    end % END IF excluding bad newest fields

    %     mydate(ifile) = data.header.date_jd;
    %     mytemp(ifile) = data.header.ccdtemp;
    %     mymean(ifile) = data.stats.calmean;
    %     myref(ifile) = data.ref.mean;
    %     myeng(ifile) = data.ref.engmean;
    %     %myisl(ifile) = data.isl.trimean;
    %
    %     if 0
    %     figure(1); clf
    %     imagesc(data.image.calimage)
    %     caxis([0,3000])
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	'(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_raw.png',data.header.timestamp));
    %     figure(2); clf
    %     imagesc(data.mask.mask)
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	  '(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_mask.png',data.header.timestamp));
    %     figure(3); clf
    %     imagesc(data.image.calimage.*~data.mask.mask)
    %     caxis([0,1500])
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	'(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_masked.png',data.header.timestamp));
    %       end


    %dbstop


% end %END IF single file

totalghosts;
totalrealghosts;
max(ghostdist);

% figure(1);
% plot(oddMinusEvenref);
% xlabel('Image Number');
% ylabel('Mean(Odd Col - Even Col) [DN]')
% title('Raw Image before Correction')
%
% figure(2);
% plot(oddMinusEvenrefFixed);
% xlabel('Image Number');
% ylabel('Mean(Odd Col - Even Col) [DN]')
% title('Raw Image after Correction')
%
% figure(3);
% plot(evenMinusOdd);
% xlabel('Image Number');
% ylabel('Mean(Even Col - Odd Col) [DN]')
% title('Cal Image before Correction')
%
% figure(4);
% plot(evenMinusOddFixed);
% xlabel('Image Number');
% ylabel('Mean(Even Col - Odd Col) [DN]')
% title('Cal Image after Correction')

%figure(4);
%plot(mydate - data.header.launch_jd,mytemp,'o');
timeCntr = timeCntr + 1;
timePer = toc;
end
