%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_fix_files.m
%%  Jun 2016, MZ
%%  This program takes in Chi's aligned astrometry information and the image
%%  data and appends them for later analysis.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nh_fix_files()

%***CAUTION: RUNNING THIS ON ESTABLISHED DATA OVERWRITES DATA MAT FILES***
% Load old data files and save if not starting from scratch

% get the paths appropriate to this system
paths = get_paths_old_ghosts();

% get a list of all the image files in the *.fit format
imagedata = dir(sprintf('%s*.fit',paths.imagedir));
% get a list of all the astrometry files in the .mat format
astrodata = dir(sprintf('%s*.mat',paths.astrodir));

% as a simple check, make sure the number of files in each is equal,
% and warn me if not
if size(imagedata) ~= size(astrodata)
    warning('Number of data images does not match number of astrometry files.')
end

% loop through all the image files
for ifile=1:size(imagedata)
    sprintf('On file number %d',ifile)
    % read in the .fits data for the ith file
    image_i = fitsread(sprintf('%s%s',...
        paths.imagedir,imagedata(ifile).name));
    
    % the astrometry file name uses only a subset of the full file name
    % specification, so strip that out
    subfile = imagedata(ifile).name(8:end);
    
    % load up the astrometry information
    load(sprintf('%s%s.mat',paths.astrodir,subfile));
    
    % If modifying established data (and not re-running from scratch), load data file here
%     load(sprintf('%s%s',datadir,datafiles(ifile).name));
    
    % append the image data to the data structure
    data.data = image_i;
    
    %create astrometry type object to use for pix2radec funtion.
    info = fitsinfo(sprintf('%s%s',...
        paths.imagedir,imagedata(ifile).name));
    astrom = fits_to_astrometry(info);
    data.astrom = astrom;
    
    % strip out the observation time stamp and add it to the header information
    timestamp = data.astrometry.id_spcutcjd(4:end);
    data.header.timestamp = timestamp;
    
    % find the field number for this observation
    fieldnum = get_field(data.astrometry.id_crval1,data.astrometry.id_crval2,paths);
    data.header.fieldnum = fieldnum;
    
    data.header.rawfile = subfile;
    
    % now add reference pixel information
    if strcmp(subfile(1:7),'jupiter')
        darkpre = subfile(9:23);
        darkpst = subfile(25:45);
    elseif strcmp(subfile(1:6),'cruise') || strcmp(subfile(1:6),'launch')
        darkpre = subfile(8:22);
        darkpst = subfile(24:44);
    elseif strcmp(subfile(1:5),'pluto')
        darkpre = subfile(7:21);
        darkpst = subfile(23:43);
    end
    
    darkstring = sprintf('%s%s/%seng_1.fit',paths.engdir,darkpre,darkpst);
%     darkstring = sprintf('%s%s/%seng.fit',paths.engdir,darkpre,darkpst);
    
    
    data.ref.file = darkstring;
    
    ref_i = fitsread(sprintf('%s',darkstring));
    
    data.ref.eng = ref_i;
    data.ref.line = ref_i(:,257);
    whpl = ~isnan(data.ref.line);
    data.ref.mean = mean(ref_i(whpl,257));
    data.ref.std = std(ref_i(whpl,257));
    
    % save the resulting data file
    save(sprintf('%s%s.mat',paths.datadir,timestamp),'data');
    
end

% and we're out

%end
