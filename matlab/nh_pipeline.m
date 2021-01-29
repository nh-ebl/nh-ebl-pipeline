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
paths = get_paths_new();

%% Set user parameters here

%set portion of pipeline to run
procstepflag = 1;

%use USNOB1 or Gaia for star masking
use_gaia = 1; %0 for USNOB1
new_star_mask = 1; %0 for skipping star masking and using existing mask

%use Gaia or just Trilegal for ISL calc
tri_gaia = 0; %0 for only Trilegal, 1 for Gaia and Trilegal

%keep or omit galaxies from Gaia data
gals = 0; %0 for galaxies removed, 1 for galaxies included

%set max mag for masking (mask all stars brighter than this), also
%calculate gaia isl between this and 20
max_mag = 21 ; %17.75 old value

%set min mag for trilegal isl (calculate isl for stars fainter than this)
tri_mag = 20;

%set max mag for psf wing isl calc (from USNOB1 or Gaia based on masking catalog)
wing_mag = max_mag;

%optionally save lists of stars and flux
save_file = 0; %1 to save

%%
%load directory of saved mat files for all images
datafiles = dir(sprintf('%s*.mat',paths.datadir));

mydate = zeros(size(datafiles));
mytemp = zeros(size(datafiles));
mymean = zeros(size(datafiles));
myref = zeros(size(datafiles));
myeng = zeros(size(datafiles));
myisl = zeros(size(datafiles));

photcurr = zeros(size(datafiles));
checkmax = zeros(size(datafiles));
checkmin = zeros(size(datafiles));

% totalghosts = 0;
if procstepflag == 3
    evenMinusOdd = zeros( size(datafiles,1), 1); %preallocate
    evenMinusOddFixed = zeros( size(datafiles,1), 1); %preallocate
end

%loop through all images
for ifile=1:size(datafiles,1)
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    %load individual mat file for image
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    %print field number
    %     disp(data.header.fieldnum)
    
    if procstepflag == 1
        fprintf('Masking and ghost analysis');
        %% Masking and optical ghosts
        %save ra, dec, pix, and mag of stars in image that could be causing ghosts
%         [data, ghostcount] = nh_findghoststar(data,paths,use_gaia);
        %         totalghosts = totalghosts + ghostcount;
        %
        %manually mask portions of image
        %         data = nh_make_manmask(data,paths);
        
        %mask out stars and ghosts in image (also stat and clip mask)
        %may need to redo catalog depending on gals
        %                 catalog_data_gaia(gals)
        data = nh_makemask(data,paths,3,use_gaia,new_star_mask, max_mag, save_file);
        
        %         maskim = data.data.*~data.mask.onemask;
        %         maskim(maskim==0) = NaN;
        %         checkmax(ifile) = nanmax(nanmax(maskim));
        %         checkmin(ifile) = nanmin(nanmin(maskim));
        
        %print total masked surface brightness (sum of all masked pixels)
        %         data.cal.sbconv .* sum(sum(data.data(data.mask.onemask)./ data.astrom.exptime))
        
        %print mask fraction
        data.mask.maskfrac
        
        %overwrite data file with changes
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 2
        fprintf('Meta data and dark current');
        %% Meta data and dark current
        %write header values, constants, coordinates, and dark reference pixel data to data file
        %(uses mask, may need to redo mask first)
        %changes values used in calibration
        data = nh_add_meta(data);
        
        %overwrite data file with changes
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 3
        fprintf('Jail bars and calibration');
        %% Jail bar correction and calibration
        %perform jail bar correction on data.data - need mask first, but can redo from original data.data
        %         [data, evenMinusOdd(ifile), evenMinusOddFixed(ifile)] = nh_jail_bars(data);
        
        %save calibrated image in surface brightness units (nh_calcref saves refcorr which is used here)
        %must redo calibration after jail bar correction, but nh_calcref must be run first
        data = nh_calibrate(data);
        
        %overwrite data file with changes
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 4
        fprintf('ISL calculation');
        %% ISL calculation
        
        %calculate ISL from USNOB1 and Trilegal
        %may need to redo catalog depending on gals
        catalog_data_gaia(gals)
        %         data = nh_calcisl(data, paths, use_gaia, tri_gaia, tri_mag, wing_mag, max_mag, save_file);
        
        %         wing = data.isl.usnowing
        %         triisl = data.isl.trimean
        %         gaiaisl = data.isl.gaiamean
        
        %overwrite data file with changes
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 5
        fprintf('DGL calculation');
        %% DGL calculation
        %         if data.header.fieldnum == 8
        %             fprintf('its here');
        %         end
        
        %calculate DGL using Planck (or IRIS) data
        data = nh_calcdgl(data, paths);
        
        %overwrite data file with changes
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
        
    end
    
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
    
    
end

% totalghosts



%figure(4);
%plot(mydate - data.header.launch_jd,mytemp,'o');


end
