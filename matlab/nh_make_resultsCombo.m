% function nh_make_resultsCombo()
close all
clear variables

% make sure parpool is the right type
parpoolobj = gcp('nocreate'); % check for thread pool, which can't use the load call
if isa(parpoolobj,'parallel.ThreadPool')
    delete(parpoolobj); %threads can't use load and will error out
end

% Set all field colors
PLOT_color = {'#1F77B4','#FF7F0E','#7256C1','#2CA02C','#D81B60','#656565',...
              '#0099A9','#3E0402','#480094','#FF4545','#004D40','#9A1A00',...
              '#224EB3','#249A7D','#FF45CA','#38FF32','#AF5D03','#00004E','#5F0038'};

% Load DGL params 
dglparams = nh_get_dgl_params();

% Set user flags
want_errmags = 0;
want_errpsf = 0;
sub_tri = 1; % 1 = subtract trilegal/psf wings, 0 = no trilegal/psf wings sub
FLG_saveGoodDataFiles = 1; %saves all used data files ("good" ones) in a separate directory
FLG_saveGoodDataFiles_dir = '/data/symons/nh_pds_data'; %directory to save "good" files to

% Import paths for data location
paths = get_paths_old();
npaths = get_paths_new();
lpaths = get_paths_lauer();
wpaths = get_paths_newest();

% Set exclude time
excl_time = 150; %seconds from start of sequence to exclude

%old light data
lightfiles = dir(sprintf('%s*.mat',paths.datadir));
%new light data
nlightfiles = dir(sprintf('%s*.mat',npaths.datadir));
%lauer light data
llightfiles = dir(sprintf('%s*.mat',lpaths.datadir));
%newest light data
wlightfiles = dir(sprintf('%s*.mat',wpaths.datadir));

%Determine which light files correspond to 'good' fields
isgoodold = zeros(numel(lightfiles),1);
isgoodnew = zeros(numel(nlightfiles),1);
isgoodlauer = zeros(numel(llightfiles),1);
isgoodnewest = zeros(numel(wlightfiles),1);

%These field numbers were previously determined by going through all the
%files
oldgoodfields = [5,3,6,7];
% oldgoodfields = [3,6,7];
% oldgoodfields = []; % Set good fields to none to skip data set
newgoodfields = [8,5,6,7];
% newgoodfields = [8,5,6,7];
newgood_exlude_enable = true; %enables skipping of sequences
% newgoodfields = [];
lauergoodfields = [1,2,3,4,5,6,7];
% lauergoodfields = [2];
lauer_exlude_enable = true; %enables skipping of new sequences
% lauergoodfields = [];
newestgoodfields = flip([2,4,6,7]); %[2,4,6,7]; %[2,4,5,6,7,12,15,16,17,19,20,23]; - old set before exclusion, use for cam cut test
% newestgoodfields = [6,2];
% newestgoodfields = flip([2,4,5,6,7,12,15,16,17,19,20,23]);
newest_exlude_enable = true; %enables skipping of new sequences
% newestgoodfields = [];


% Calculate total fields
goodfiles = [oldgoodfields,newgoodfields,newestgoodfields,lauergoodfields];
zones = cumsum([1,length(oldgoodfields),length(newgoodfields),length(newestgoodfields),length(lauergoodfields)]); %Record the seperate zones and their orders
unique_fields = [1:length(goodfiles)];

num_fields = length(goodfiles);
field_snaps_data = zeros(256,256,num_fields); %records the 1st good image for each unique field
field_snaps_cal = zeros(256,256,num_fields); %records the 1st good image for each unique field
field_snaps_mask = true(num_fields,1); %to record which field snaps have been recorded
get3DInd = @(tf,pg)find(tf) + numel(tf)*(find(pg)-1); %makes 2d logical arrays work in 3d array

%Check for old light files
zones_indexer = 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
zones_snaps_mask = field_snaps_mask(zones(zones_indexer):zones(zones_indexer+1)-1); %current zone's snaps mask
zones_snaps_data = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
zones_snaps_cal = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
for ifile=1:numel(lightfiles)
    datatemp = load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if any(data.header.fieldnum == oldgoodfields)
        isgoodold(ifile) = 1;
    end
    if isgoodold(ifile) && any(data.header.fieldnum == zones_goodfiles(zones_snaps_mask))
        zones_snaps_data(:,:,data.header.fieldnum == zones_goodfiles) = data.data;
        zones_snaps_cal(:,:,data.header.fieldnum == zones_goodfiles) = data.image.calimage;
        zones_snaps_cal(get3DInd(data.mask.onemask,data.header.fieldnum == zones_goodfiles)) = nan; %two step to set masked values as nan
        zones_snaps_mask(data.header.fieldnum == zones_goodfiles) = 0; %turn off, got it
    end
end
field_snaps_data(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_data; %copy it in
field_snaps_cal(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_cal;
oldgood_exclude_logical = false(size(isgoodold)); %not used
old_goodfilesNum = sum(isgoodold&~oldgood_exclude_logical);
disp(['Old|Total files used after isgoodold (good field check, ',num2str(sum(isgoodold)),') and ~oldgood_exclude (time skip, NOT USED): ',num2str(old_goodfilesNum),'/',num2str(numel(lightfiles))])

%Check for new light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = excl_time; %s, time to skip at start of sequence
newgood_exclude = zeros(numel(nlightfiles),1); %will fill up with new sequences to ignore
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
zones_snaps_mask = field_snaps_mask(zones(zones_indexer):zones(zones_indexer+1)-1); %current zone's snaps mask
zones_snaps_data = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
zones_snaps_cal = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
for ifile=1:numel(nlightfiles)
    datatemp = load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == newgoodfields)
        isgoodnew(ifile) = 1;
    end
    if newgood_exlude_enable
        if( ifile == 1 )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to skip dynamically
        end
        %skip a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            newgood_exclude(ifile) = 1; %set this to exlude
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            newgood_exclude(ifile) = 1; %set this to exlude
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to skip
        end
    end
    if isgoodnew(ifile) && ~newgood_exclude(ifile) && any(data.header.fieldnum == zones_goodfiles(zones_snaps_mask))
        zones_snaps_data(:,:,data.header.fieldnum == zones_goodfiles) = data.data;
        zones_snaps_cal(:,:,data.header.fieldnum == zones_goodfiles) = data.image.calimage;
        zones_snaps_cal(get3DInd(data.mask.onemask,data.header.fieldnum == zones_goodfiles)) = nan; %two step to set masked values as nan
        zones_snaps_mask(data.header.fieldnum == zones_goodfiles) = 0; %turn off, got it
    end
end
field_snaps_data(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_data; %copy it in
field_snaps_cal(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_cal;
newgood_exclude_logical = newgood_exclude;
newgood_exclude = find(newgood_exclude); %get it into indexes
new_goodfilesNum = sum(isgoodnew&~newgood_exclude_logical);
disp(['New|Total files used after isgoodnew (good field check, ',num2str(sum(isgoodnew)),') and ~newgood_exclude (time skip, ',num2str(sum(~newgood_exclude_logical)),'): ',num2str(new_goodfilesNum),'/',num2str(numel(nlightfiles)),' [there may be overlap in the two values]'])

%Check for newest light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = excl_time; %s, time to skip at start of sequence
newest_exclude = zeros(numel(wlightfiles),1); %will fill up with new sequences to ignore
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
zones_snaps_mask = field_snaps_mask(zones(zones_indexer):zones(zones_indexer+1)-1); %current zone's snaps mask
zones_snaps_data = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
zones_snaps_cal = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
for ifile=1:numel(wlightfiles)
    load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
    if sum(data.header.fieldnum == newestgoodfields)
        isgoodnewest(ifile) = 1;
    end
    if newest_exlude_enable
        if( ifile == 1 )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to skip dynamically
        end
        %skip a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            newest_exclude(ifile) = 1; %set this to exlude
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
%             disp(['Current req ID: ',data.astrom.reqid,' Field #: ',num2str(data.header.fieldnum),' Excluded so far: ',num2str(sum(newest_exclude))])
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            newest_exclude(ifile) = 1; %set this to exlude
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to skip
        end
    end
    if isgoodnewest(ifile) && ~newest_exclude(ifile) && any(data.header.fieldnum == zones_goodfiles(zones_snaps_mask))
        zones_snaps_data(:,:,data.header.fieldnum == zones_goodfiles) = data.data;
        zones_snaps_cal(:,:,data.header.fieldnum == zones_goodfiles) = data.image.calimage;
        zones_snaps_cal(get3DInd(data.mask.onemask,data.header.fieldnum == zones_goodfiles)) = nan; %two step to set masked values as nan
        zones_snaps_mask(data.header.fieldnum == zones_goodfiles) = 0; %turn off, got it
    end
end
field_snaps_data(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_data; %copy it in
field_snaps_cal(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_cal;
newest_exclude_logical = newest_exclude;
newest_exclude = find(newest_exclude); %get it into indexes
% disp(['Total excluded: ',num2str(sum(newest_exclude))])
newest_goodfilesNum = sum(isgoodnewest&~newest_exclude_logical);
disp(['Newest|Total files used after isgoodnewest (good field check, ',num2str(sum(isgoodnewest)),') and ~newest_exclude (time skip, ',num2str(sum(~newest_exclude_logical)),'): ',num2str(newest_goodfilesNum),'/',num2str(numel(wlightfiles)),' [there may be overlap in the two values]'])

%Check for lauer light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = excl_time; %s, time to skip at start of sequence
lauer_exclude = zeros(numel(llightfiles),1); %will fill up with new sequences to ignore
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
zones_snaps_mask = field_snaps_mask(zones(zones_indexer):zones(zones_indexer+1)-1); %current zone's snaps mask
zones_snaps_data = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
zones_snaps_cal = zeros(256,256,length(zones_snaps_mask)); %records the 1st good image for each unique field
for ifile=1:numel(llightfiles)
    load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
    if sum(data.header.fieldnum == lauergoodfields)
        isgoodlauer(ifile) = 1;
    end
    if lauer_exlude_enable
        if( ifile == 1 )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to skip dynamically
        end
        %skip a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            lauer_exclude(ifile) = 1; %set this to exlude
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            lauer_exclude(ifile) = 1; %set this to exlude
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to skip
        end
    end
    if isgoodlauer(ifile) && ~lauer_exclude(ifile) && any(data.header.fieldnum == zones_goodfiles(zones_snaps_mask))
        zones_snaps_data(:,:,data.header.fieldnum == zones_goodfiles) = data.data;
        zones_snaps_cal(:,:,data.header.fieldnum == zones_goodfiles) = data.image.calimage;
        zones_snaps_cal(get3DInd(data.mask.onemask,data.header.fieldnum == zones_goodfiles)) = nan; %two step to set masked values as nan
        zones_snaps_mask(data.header.fieldnum == zones_goodfiles) = 0; %turn off, got it
    end
end
field_snaps_data(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_data; %copy it in
field_snaps_cal(:,:,zones(zones_indexer):zones(zones_indexer+1)-1) = zones_snaps_cal;
lauer_exclude_logical = lauer_exclude;
lauer_exclude = find(lauer_exclude); %get it into indexes
lauer_goodfilesNum = sum(isgoodlauer&~lauer_exclude_logical);
disp(['Lauer|Total files used after isgoodlauer (good field check, ',num2str(sum(isgoodlauer)),') and ~lauer_exclude (time skip, ',num2str(sum(~lauer_exclude_logical)),'): ',num2str(lauer_goodfilesNum),'/',num2str(numel(llightfiles)),' [there may be overlap in the two values]'])


load('run_params.mat','params')

% If method file exists, read saved text file for data_type and see which method last used
flag_method = 'new';

%Number of old and new light files corresponding to good fields
numoldlightfiles = numel(lightfiles);
numnewlightfiles = numel(nlightfiles);
numlauerlightfiles = numel(llightfiles);
numnewestlightfiles = numel(wlightfiles);
% %Number of old and new light files corresponding to good fields
% numoldlightfiles = sum(isgoodold);
% numnewlightfiles = sum(isgoodnew);
% numlauerlightfiles = sum(isgoodlauer);
% numnewestlightfiles = sum(isgoodnewest);

%Preallocate
file_num_perSet = [1:numoldlightfiles,1:numnewlightfiles,1:numnewestlightfiles,1:numlauerlightfiles]'; %the number of the file in the specific set
isgood = false((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mydate = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myfieldnum = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
unique_field_id = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mydataset = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mytarget = cellstr('');
mysun = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mydt = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myelong = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles)+numnewestlightfiles,1);
myexp = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myref = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myeng = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myohm = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
% myav = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mysig = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mysig_magerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mymasked = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mymasked_magerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
biaslevl = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
biasmthd = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
biasdiff = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
rsolar = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myb = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mytri = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mytrierr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mypsfwing = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mypsfwing_psferr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myghost = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_tot = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_toterr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_gaia = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_gaiaerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_gaia_psferrpos= zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_gaia_fluxerrpos = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_masana = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_masanaerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_masana_psferrpos = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myscattering_masana_fluxerrpos = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myisl = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mydgl = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mydglerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mycrr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mymaskmean = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myext = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myghostdiff = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myghostdifferrpos = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myghostdifferrneg = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
myghostdiffcomp = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_mydgl_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_mydglerr_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
% pltr_myohmim_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles),1);
pltr_mydgl_iris = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_mydglerr_iris = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
% pltr_myohmim_iris = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles),1);
pltr_mydgl_iris_sfd = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_mydglerr_iris_sfd = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100m_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_beta_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_temp_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_tau_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_combo_planck = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100m_iris = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_iris = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100m_iris_sfd = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_my100merr_iris_sfd = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_mydgl_nh = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
pltr_mydglerr_nh = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mygalerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mygal_mean = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mymagerr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
mypsferr = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
maskfrac = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);

%For old data files
fprintf('Loading old data \n');
jfile_offset = 0; %offset to keep ifile in line with other file sets
% k = 1;
zones_indexer = 1; %Keeps count of the current zone
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
parfor jfile=1+jfile_offset:numel(lightfiles)+jfile_offset
    %If file is for a good field, load and save values
    ifile = jfile-jfile_offset; %lets parfor work by letting jfile iterate instead of ifile
    if isgoodold(ifile) == 1
        %Load data files
        dataz = load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(lightfiles,1)));
        data = dataz.data; %get the data out, makes parfor work
        
        if( FLG_saveGoodDataFiles == 1 )
            [copy_status, copy_msg] = copyfile(sprintf('%s%s',paths.datadir,lightfiles(ifile).name), [FLG_saveGoodDataFiles_dir,'/',lightfiles(ifile).name]);
            if( copy_status == 0 )
                disp(['WARNING: file ',sprintf('%s%s',paths.datadir,lightfiles(ifile).name),' failed to be copied to ',[FLG_saveGoodDataFiles_dir,'/',lightfiles(ifile).name],' failure message follows:'])
                disp(copy_msg)
            end
        end

        %Read in data
        mydate(jfile) = data.header.date_jd-data.header.launch_jd;

        myfieldnum(jfile) = data.header.fieldnum;
        %Assign unique field id based on field number
        unique_field_id(jfile) = find(myfieldnum(jfile) == zones_goodfiles) + zones(zones_indexer)-1; %Counts from 1 to length(goodfiles), like unique_fields
        %         if jfile == 1
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) ~= myfieldnum(jfile-1)
        %             k = k + 1;
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) == myfieldnum(jfile-1)
        %             unique_field_id(jfile) = unique_fields(k);
        %         end
        mydataset(jfile) = zones_indexer; %record the dataset of the data (via the zone, 1 == old, 2 == new, 3 == lauer)
        
        mytarget{jfile} = data.header.date_cal;
        
        mysun(jfile) = data.distance.d_sun_NH;

        mydt(jfile) = data.distance.d_NH_target;%./1.496e8;;

        myelong(jfile) = data.coords.sol_elon;

        myb(jfile) = data.coords.galactic(2);

        myexp(jfile) = round(data.header.exptime);

        myref(jfile) = data.ref.mean;

        myeng(jfile) = data.ref.engmean;

        %         myohm(jfile) = data.dgl.ohmmean;

        %         myav(jfile) = 0.4.*data.header.A_V;

        myext(jfile) = data.ext.mean_flux_rat;

        % Reference-corrected mean of masked image
        if want_errmags == 1
            mysig(jfile) = data.err_mags.stats.corrmean;
        else
            mysig(jfile) = data.stats.corrmean;
        end
        mysig_magerr(jfile) = data.err_mags.stats.corrmean;

        mymaskmean(jfile) = data.stats.maskmean;

        maskfrac(jfile) = data.mask.maskfrac;

        % Bias level from header
        %     biaslevl(ifile) = data.astrom.biaslevl;
        % Bias method from header (mean = 1, median = 2) and difference to
        % actual mean/median
        %     if strcmp(data.astrom.biasmthd,'mean of dark column data') == 1
        %         biasmthd(ifile) = 1;
        %         biasdiff(ifile) = data.astrom.biaslevl - mean(data.ref.line);
        %     elseif strcmp(data.astrom.biasmthd,'median of dark column data') == 1
        %         biasmthd(ifile) = 2;
        %         biasdiff(ifile) = data.astrom.biaslevl - median(data.ref.line);
        %     else
        %         fprintf('Absolute panic: bias method not recognized!')
        %     end
        % RSOLAR from header
        rsolar(jfile) = data.astrom.rsolar;

        % Reference-corrected most probable value of masked image
        %     mysig(ifile) = data.stats.corrmostprob;

        % Sum of image brightness of *masked* stars - how bright is the
        % masked-out portion of each image per pixel
        if want_errmags == 1
            starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        else
            starmaskfrac = sum(data.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.image.calimage.*data.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        end
        starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
        mymasked_magerr(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);

        mytri(jfile) = data.isl.trimean;
        mytrierr(jfile) = data.isl.trierr;

        if want_errpsf == 1
            mypsfwing(jfile) = data.err_psf.isl.usnowing;
        else
            mypsfwing(jfile) = data.isl.usnowing;
        end
        mypsfwing_psferr(jfile) = data.err_psf.isl.usnowing;
        
        if sub_tri == 1
            myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);
        elseif sub_tri == 0
            myisl(jfile) = 0;
        end

        pltr_mydgl_planck(jfile) = data.dgl.dglmean_planck;
        pltr_mydglerr_planck(jfile) = data.dgl.dglerr_planck;
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_mydgl_iris(jfile) = data.dgl.dglmean_iris;
        pltr_mydglerr_iris(jfile) = data.dgl.dglerr_iris;
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_mydgl_iris_sfd(jfile) = data.dgl.dglmean_iris_sfd;
        pltr_mydglerr_iris_sfd(jfile) = data.dgl.dglerr_iris_sfd;
        pltr_my100m_planck(jfile) = data.dgl.onehundomean_planck;
        pltr_my100merr_planck(jfile) = data.dgl.ohmstd_planck;
        pltr_my100merr_beta_planck(jfile) = data.dgl.sem_planck_mc_beta;
        pltr_my100merr_temp_planck(jfile) = data.dgl.sem_planck_mc_temp;
        pltr_my100merr_tau_planck(jfile) = data.dgl.sem_planck_mc_tau;
        % Combine temp and beta errors in quad because they're usually equal in magnitude
        pltr_my100merr_combo_planck(jfile) = sqrt(data.dgl.sem_planck_mc_temp^2 + data.dgl.sem_planck_mc_beta^2);
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_my100m_iris(jfile) = data.dgl.onehundomean_iris;
        pltr_my100merr_iris(jfile) = std(std(data.dgl.ohmim_iris-0.48));
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_my100m_iris_sfd(jfile) = data.dgl.onehundomean_iris_sfd;
        pltr_my100merr_iris_sfd(jfile) = data.dgl.ohmstd_iris_sfd;
        pltr_mydgl_nh(jfile) = data.dgl.ohmmean_nh;
        pltr_mydglerr_nh(jfile) = data.dgl.ohmstd_nh;

        % Calculate mean difference between real cal mean and err sims
        mygalerr(jfile) = mean(data.stats.calmean - data.err_gals.stats.calmean_mc(:));
        mygal_mean(jfile) = mean(data.err_gals.stats.calmean_mc(:));
        mymagerr(jfile) = data.stats.calmean - data.err_mags.stats.calmean;
        mypsferr(jfile) = data.isl.usnowing - data.err_psf.isl.usnowing;

        if want_errmags == 1
            myerr(jfile) = data.err_mags.stats.correrr;
        else
            myerr(jfile) = data.stats.correrr;
        end

        if want_errmags == 1
            mycrr(jfile) = data.err_mags.ref.bias;
        else
            mycrr(jfile) = data.ref.bias;
        end

        if strcmp(flag_method,'new') == 1
            if want_errmags == 1
                myghostdiff(jfile) = data.err_mags.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.err_mags.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.err_mags.ghost.diffusesuberrneg;
            else
                myghostdiff(jfile) = data.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.ghost.diffusesuberrneg;
                %             myghostdiffcomp(jfile) = data.ghost.diffusediff;
                myscattering_tot(jfile) = data.scattered.gaiasum + data.scattered.masanasum;
                myscattering_gaia(jfile) = data.scattered.gaiasum;
                myscattering_gaia_psferrpos(jfile) = data.scattered.gaiasum_psferrmax - data.scattered.gaiasum;
                myscattering_gaia_fluxerrpos(jfile) = data.scattered.gaiasum_fluxerrmax - data.scattered.gaiasum;
                myscattering_gaiaerr(jfile) = sqrt(myscattering_gaia_psferrpos(jfile)^2 + myscattering_gaia_fluxerrpos(jfile)^2);
                myscattering_masana(jfile) = data.scattered.masanasum;
                myscattering_masana_psferrpos(jfile) = data.scattered.masanasum_psferrmax - data.scattered.masanasum;
                myscattering_masana_fluxerrpos(jfile) = data.scattered.masanasum_fluxerrmax - data.scattered.masanasum;
                myscattering_masanaerr(jfile) = sqrt(myscattering_masana_psferrpos(jfile)^2 + myscattering_masana_fluxerrpos(jfile)^2);
                myscattering_toterr(jfile) = myscattering_gaiaerr(jfile) + myscattering_masanaerr(jfile);
            end
        end


        % If struct field data.header.bad exists, check if file is good or
        % bad - only mark as good if not bad
        if isfield(data.header,'bad')
            % Only use good images
            if data.header.bad == 0
                isgood(jfile) = 1;
                mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                    mydate(jfile),myfieldnum(jfile),...
                    mytarget{jfile},mysun(jfile),myelong(jfile),...
                    mycrr(jfile),...
                    mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                    mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

                disp(mystring)

                image = data.image.calimage;
                % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
                parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
                parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
            else
                if( FLG_saveGoodDataFiles == 1 )
                    delete([FLG_saveGoodDataFiles_dir,'/',lightfiles(ifile).name]) % Remove because it actually was a bad file all along
                end
            end
            % If struct field does not exist, files have not been marked good
            % or bad - assume all good
        else
            isgood(jfile) = 1;
            mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                mydate(jfile),myfieldnum(jfile),...
                mytarget{jfile},mysun(jfile),myelong(jfile),...
                mycrr(jfile),...
                mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

            disp(mystring)

            image = data.image.calimage;
            % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
            parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
            parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
        end
    end
end
disp(['Total files being used for old data: ',num2str(sum(isgood(1+jfile_offset:numel(lightfiles)+jfile_offset) & isgoodold & ~oldgood_exclude_logical)),'/',num2str(numel(lightfiles))])
disp(['isgood (header.bad [also incl. field check and time skip]): ',num2str(sum(isgood(1+jfile_offset:numel(lightfiles)+jfile_offset))),' | ',...
    'isgoodold (good field check): ',num2str(sum(isgoodold)),' | ',...
    'oldgood_exclude (time skip [NOT USED]): ',num2str(sum(~oldgood_exclude_logical))])

%For new data files
jfile_offset = numoldlightfiles; %offset to keep ifile in line with other file sets
fprintf('Loading new data \n');
badcntr = 0;
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
parfor jfile=1+jfile_offset:numel(nlightfiles)+jfile_offset
    %If file is for a good field, load and save values
    ifile = jfile-jfile_offset; %lets parfor work by letting jfile iterate instead of ifile
    if (isgoodnew(ifile) == 1) && ~any(newgood_exclude == ifile)
        %Load data files
        dataz = load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(nlightfiles,1)));
        data = dataz.data; %get the data out, makes parfor work
        
        if( FLG_saveGoodDataFiles == 1 )
            [copy_status, copy_msg] = copyfile(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name), [FLG_saveGoodDataFiles_dir,'/',nlightfiles(ifile).name]);
            if( copy_status == 0 )
                disp(['WARNING: file ',sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name),' failed to be copied to ',[FLG_saveGoodDataFiles_dir,'/',nlightfiles(ifile).name],' failure message follows:'])
                disp(copy_msg)
            end
        end

        %Read in data
        mydate(jfile) = data.header.date_jd-data.header.launch_jd;

        myfieldnum(jfile) = data.header.fieldnum;
        %Assign unique field id based on field number
        unique_field_id(jfile) = find(myfieldnum(jfile) == zones_goodfiles) + zones(zones_indexer)-1; %Counts from 1 to length(goodfiles), like unique_fields
        %         if jfile == 1
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) ~= myfieldnum(jfile-1)
        %             k = k + 1;
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) == myfieldnum(jfile-1)
        %             unique_field_id(jfile) = unique_fields(k);
        %         end
        mydataset(jfile) = zones_indexer; %record the dataset of the data (via the zone, 1 == old, 2 == new, 3 == lauer)

        mytarget{jfile} = data.header.date_cal;

        mysun(jfile) = data.distance.d_sun_NH;

        mydt(jfile) = data.distance.d_NH_target;%./1.496e8;;

        myelong(jfile) = data.coords.sol_elon;

        myb(jfile) = data.coords.galactic(2);

        myexp(jfile) = round(data.header.exptime);

        myref(jfile) = data.ref.mean;

        myeng(jfile) = data.ref.engmean;

        maskfrac(jfile) = data.mask.maskfrac;

        %         myohm(jfile) = data.dgl.ohmmean;

        %         myav(jfile) = 0.4.*data.header.A_V;

        myext(jfile) = data.ext.mean_flux_rat;

        % Reference-corrected mean of masked image
        if want_errmags == 1
            mysig(jfile) = data.err_mags.stats.corrmean;
        else
            mysig(jfile) = data.stats.corrmean;
        end
        mysig_magerr(jfile) = data.err_mags.stats.corrmean;

        mymaskmean(jfile) = data.stats.maskmean;

        % Bias level from header
        %     biaslevl(ifile) = data.astrom.biaslevl;
        % Bias method from header (mean = 1, median = 2) and difference to
        % actual mean/median
        %     if strcmp(data.astrom.biasmthd,'mean of dark column data') == 1
        %         biasmthd(ifile) = 1;
        %         biasdiff(ifile) = data.astrom.biaslevl - mean(data.ref.line);
        %     elseif strcmp(data.astrom.biasmthd,'median of dark column data') == 1
        %         biasmthd(ifile) = 2;
        %         biasdiff(ifile) = data.astrom.biaslevl - median(data.ref.line);
        %     else
        %         fprintf('Absolute panic: bias method not recognized!')
        %     end
        % RSOLAR from header
        rsolar(jfile) = data.astrom.rsolar;

        % Reference-corrected most probable value of masked image
        %     mysig(ifile) = data.stats.corrmostprob;

        % Sum of image brightness of *masked* stars - how bright is the
        % masked-out portion of each image per pixel
        if want_errmags == 1
            starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        else
            starmaskfrac = sum(data.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.image.calimage.*data.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        end
        starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
        mymasked_magerr(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);

        mytri(jfile) = data.isl.trimean;
        mytrierr(jfile) = data.isl.trierr;

        if want_errpsf == 1
            mypsfwing(jfile) = data.err_psf.isl.usnowing;
        else
            mypsfwing(jfile) = data.isl.usnowing;
        end
        mypsfwing_psferr(jfile) = data.err_psf.isl.usnowing;

        if sub_tri == 1
            myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);
        elseif sub_tri == 0
            myisl(jfile) = 0;
        end

        pltr_mydgl_planck(jfile) = data.dgl.dglmean_planck;
        pltr_mydglerr_planck(jfile) = data.dgl.dglerr_planck;
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_mydgl_iris(jfile) = data.dgl.dglmean_iris;
        pltr_mydglerr_iris(jfile) = data.dgl.dglerr_iris;
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_mydgl_iris_sfd(jfile) = data.dgl.dglmean_iris_sfd;
        pltr_mydglerr_iris_sfd(jfile) = data.dgl.dglerr_iris_sfd;
        pltr_my100m_planck(jfile) = data.dgl.onehundomean_planck;
        pltr_my100merr_planck(jfile) = data.dgl.ohmstd_planck;
        pltr_my100merr_beta_planck(jfile) = data.dgl.sem_planck_mc_beta;
        pltr_my100merr_temp_planck(jfile) = data.dgl.sem_planck_mc_temp;
        pltr_my100merr_tau_planck(jfile) = data.dgl.sem_planck_mc_tau;
        % Combine temp and beta errors in quad because they're usually equal in magnitude
        pltr_my100merr_combo_planck(jfile) = sqrt(data.dgl.sem_planck_mc_temp^2 + data.dgl.sem_planck_mc_beta^2);
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_my100m_iris(jfile) = data.dgl.onehundomean_iris;
        pltr_my100merr_iris(jfile) = std(std(data.dgl.ohmim_iris-0.48));
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_my100m_iris_sfd(jfile) = data.dgl.onehundomean_iris_sfd;
        pltr_my100merr_iris_sfd(jfile) = data.dgl.ohmstd_iris_sfd;
        pltr_mydgl_nh(jfile) = data.dgl.ohmmean_nh;
        pltr_mydglerr_nh(jfile) = data.dgl.ohmstd_nh;

        % Calculate mean difference between real cal mean and err sims
        mygalerr(jfile) = mean(data.stats.calmean - data.err_gals.stats.calmean_mc(:));
        mygal_mean(jfile) = mean(data.err_gals.stats.calmean_mc(:));
        mymagerr(jfile) = data.stats.calmean - data.err_mags.stats.calmean;
        mypsferr(jfile) = data.isl.usnowing - data.err_psf.isl.usnowing;

        if want_errmags == 1
            myerr(jfile) = data.err_mags.stats.correrr;
        else
            myerr(jfile) = data.stats.correrr;
        end

        if want_errmags == 1
            mycrr(jfile) = data.err_mags.ref.bias;
        else
            mycrr(jfile) = data.ref.bias;
        end

        if strcmp(flag_method,'new') == 1
            if want_errmags == 1
                myghostdiff(jfile) = data.err_mags.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.err_mags.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.err_mags.ghost.diffusesuberrneg;
            else
                myghostdiff(jfile) = data.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.ghost.diffusesuberrneg;
                %             myghostdiffcomp(jfile) = data.ghost.diffusediff;
                myscattering_tot(jfile) = data.scattered.gaiasum + data.scattered.masanasum;
                myscattering_gaia(jfile) = data.scattered.gaiasum;
                myscattering_gaia_psferrpos(jfile) = data.scattered.gaiasum_psferrmax - data.scattered.gaiasum;
                myscattering_gaia_fluxerrpos(jfile) = data.scattered.gaiasum_fluxerrmax - data.scattered.gaiasum;
                myscattering_gaiaerr(jfile) = sqrt(myscattering_gaia_psferrpos(jfile)^2 + myscattering_gaia_fluxerrpos(jfile)^2);
                myscattering_masana(jfile) = data.scattered.masanasum;
                myscattering_masana_psferrpos(jfile) = data.scattered.masanasum_psferrmax - data.scattered.masanasum;
                myscattering_masana_fluxerrpos(jfile) = data.scattered.masanasum_fluxerrmax - data.scattered.masanasum;
                myscattering_masanaerr(jfile) = sqrt(myscattering_masana_psferrpos(jfile)^2 + myscattering_masana_fluxerrpos(jfile)^2);
                myscattering_toterr(jfile) = myscattering_gaiaerr(jfile) + myscattering_masanaerr(jfile);
            end
        end


        % If struct field data.header.bad exists, check if file is good or
        % bad - only mark as good if not bad
        if isfield(data.header,'bad')
            % Only use good images
            if data.header.bad == 0
                isgood(jfile) = 1;
                mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                    mydate(jfile),myfieldnum(jfile),...
                    mytarget{jfile},mysun(jfile),myelong(jfile),...
                    mycrr(jfile),...
                    mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                    mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

                disp(mystring)

                image = data.image.calimage;
                % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
                parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
                parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
            else
                disp(['header.bad WAS BAD on file (',num2str(ifile),'/',num2str(size(nlightfiles,1)),') ',sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name)])
                badcntr = badcntr + 1;
                if( FLG_saveGoodDataFiles == 1 )
                    delete([FLG_saveGoodDataFiles_dir,'/',nlightfiles(ifile).name]) % Remove because it actually was a bad file all along
                end
            end
            % If struct field does not exist, files have not been marked good
            % or bad - assume all good
        else
            isgood(jfile) = 1;
            mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                mydate(jfile),myfieldnum(jfile),...
                mytarget{jfile},mysun(jfile),myelong(jfile),...
                mycrr(jfile),...
                mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

            disp(mystring)

            image = data.image.calimage;
            % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
            parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
            parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
        end
    end
end
disp(['Bad cntr: ',num2str(badcntr)])
disp(['Total files being used for new data: ',num2str(sum(isgood(1+jfile_offset:numel(nlightfiles)+jfile_offset) & isgoodnew & ~newgood_exclude_logical)),'/',num2str(numel(nlightfiles))])
disp(['isgood (header.bad [also incl. field check and time skip]): ',num2str(sum(isgood(1+jfile_offset:numel(nlightfiles)+jfile_offset))),' | ',...
    'isgoodnew (good field check): ',num2str(sum(isgoodnew)),' | ',...
    'newgood_exclude (time skip): ',num2str(sum(~newgood_exclude_logical))])

%For newest data files
fprintf('Loading newest data \n');
jfile_offset = numoldlightfiles+numnewlightfiles; %offset to keep ifile in line with other file sets
badcntr = 0;
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
parfor jfile=1+jfile_offset:numel(wlightfiles)+jfile_offset
    %If file is for a good field and not excluded, load and save values
    ifile = jfile-jfile_offset; %lets parfor work by letting jfile iterate instead of ifile
    if (isgoodnewest(ifile) == 1) && ~any(newest_exclude == ifile)

        %Load data files
        dataz = load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(wlightfiles,1)));
        data = dataz.data; %get the data out, makes parfor work

        if( FLG_saveGoodDataFiles == 1 )
            [copy_status, copy_msg] = copyfile(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name), [FLG_saveGoodDataFiles_dir,'/',wlightfiles(ifile).name]);
            if( copy_status == 0 )
                disp(['WARNING: file ',sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name),' failed to be copied to ',[FLG_saveGoodDataFiles_dir,'/',wlightfiles(ifile).name],' failure message follows:'])
                disp(copy_msg)
            end
        end

        %Read in data
        mydate(jfile) = data.header.date_jd-data.header.launch_jd;

        myfieldnum(jfile) = data.header.fieldnum;
        %Assign unique field id based on field number
        unique_field_id(jfile) = find(myfieldnum(jfile) == zones_goodfiles) + zones(zones_indexer)-1; %Counts from 1 to length(goodfiles), like unique_fields
        %         if jfile == 1
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) ~= myfieldnum(jfile-1)
        %             k = k + 1;
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) == myfieldnum(jfile-1)
        %             unique_field_id(jfile) = unique_fields(k);
        %         end
        mydataset(jfile) = zones_indexer; %record the dataset of the data (via the zone, 1 == old, 2 == new, 3 == lauer)

        mytarget{jfile} = data.header.date_cal;

        mysun(jfile) = data.distance.d_sun_NH;

        mydt(jfile) = data.distance.d_NH_target;%./1.496e8;;

        myelong(jfile) = data.coords.sol_elon;

        myb(jfile) = data.coords.galactic(2);

        myexp(jfile) = round(data.header.exptime);

        myref(jfile) = data.ref.mean;

        maskfrac(jfile) = data.mask.maskfrac;


        myeng(jfile) = data.ref.engmean;

        %         myohm(jfile) = data.dgl.ohmmean;

        %         myav(jfile) = 0.4.*data.header.A_V;

        myext(jfile) = data.ext.mean_flux_rat;

        % Reference-corrected mean of masked image
        if want_errmags == 1
            mysig(jfile) = data.err_mags.stats.corrmean;
        else
            mysig(jfile) = data.stats.corrmean;
        end
        mysig_magerr(jfile) = data.err_mags.stats.corrmean;

        mymaskmean(jfile) = data.stats.maskmean;

        % Bias level from header
        %     biaslevl(ifile) = data.astrom.biaslevl;
        % Bias method from header (mean = 1, median = 2) and difference to
        % actual mean/median
        %     if strcmp(data.astrom.biasmthd,'mean of dark column data') == 1
        %         biasmthd(ifile) = 1;
        %         biasdiff(ifile) = data.astrom.biaslevl - mean(data.ref.line);
        %     elseif strcmp(data.astrom.biasmthd,'median of dark column data') == 1
        %         biasmthd(ifile) = 2;
        %         biasdiff(ifile) = data.astrom.biaslevl - median(data.ref.line);
        %     else
        %         fprintf('Absolute panic: bias method not recognized!')
        %     end
        % RSOLAR from header
        rsolar(jfile) = data.astrom.rsolar;

        % Reference-corrected most probable value of masked image
        %     mysig(ifile) = data.stats.corrmostprob;

        % Sum of image brightness of *masked* stars - how bright is the
        % masked-out portion of each image per pixel
        if want_errmags == 1
            starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        else
            starmaskfrac = sum(data.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.image.calimage.*data.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        end
        starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
        mymasked_magerr(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);

        mytri(jfile) = data.isl.trimean;
        mytrierr(jfile) = data.isl.trierr;

        if want_errpsf == 1
            mypsfwing(jfile) = data.err_psf.isl.usnowing;
        else
            mypsfwing(jfile) = data.isl.usnowing;
        end
        mypsfwing_psferr(jfile) = data.err_psf.isl.usnowing;

        if sub_tri == 1
            myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);
        elseif sub_tri == 0
            myisl(jfile) = 0;
        end

        pltr_mydgl_planck(jfile) = data.dgl.dglmean_planck;
        pltr_mydglerr_planck(jfile) = data.dgl.dglerr_planck;
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_mydgl_iris(jfile) = data.dgl.dglmean_iris;
        pltr_mydglerr_iris(jfile) = data.dgl.dglerr_iris;
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_mydgl_iris_sfd(jfile) = data.dgl.dglmean_iris_sfd;
        pltr_mydglerr_iris_sfd(jfile) = data.dgl.dglerr_iris_sfd;
        pltr_my100m_planck(jfile) = data.dgl.onehundomean_planck;
        pltr_my100merr_planck(jfile) = data.dgl.ohmstd_planck;
        pltr_my100merr_beta_planck(jfile) = data.dgl.sem_planck_mc_beta;
        pltr_my100merr_temp_planck(jfile) = data.dgl.sem_planck_mc_temp;
        pltr_my100merr_tau_planck(jfile) = data.dgl.sem_planck_mc_tau;
        % Combine temp and beta errors in quad because they're usually equal in magnitude
        pltr_my100merr_combo_planck(jfile) = sqrt(data.dgl.sem_planck_mc_temp^2 + data.dgl.sem_planck_mc_beta^2);
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_my100m_iris(jfile) = data.dgl.onehundomean_iris;
        pltr_my100merr_iris(jfile) = std(std(data.dgl.ohmim_iris-0.48));
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_my100m_iris_sfd(jfile) = data.dgl.onehundomean_iris_sfd;
        pltr_my100merr_iris_sfd(jfile) = data.dgl.ohmstd_iris_sfd;
        pltr_mydgl_nh(jfile) = data.dgl.ohmmean_nh;
        pltr_mydglerr_nh(jfile) = data.dgl.ohmstd_nh;

        % Calculate mean difference between real cal mean and err sims
        mygalerr(jfile) = mean(data.stats.calmean - data.err_gals.stats.calmean_mc(:));
        mygal_mean(jfile) = mean(data.err_gals.stats.calmean_mc(:));
        mymagerr(jfile) = data.stats.calmean - data.err_mags.stats.calmean;
        mypsferr(jfile) = data.isl.usnowing - data.err_psf.isl.usnowing;

        if want_errmags == 1
            myerr(jfile) = data.err_mags.stats.correrr;
        else
            myerr(jfile) = data.stats.correrr;
        end

        if want_errmags == 1
            mycrr(jfile) = data.err_mags.ref.bias;
        else
            mycrr(jfile) = data.ref.bias;
        end

        if strcmp(flag_method,'new') == 1
            if want_errmags == 1
                myghostdiff(jfile) = data.err_mags.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.err_mags.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.err_mags.ghost.diffusesuberrneg;
            else
                myghostdiff(jfile) = data.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.ghost.diffusesuberrneg;
                %             myghostdiffcomp(jfile) = data.ghost.diffusediff;
                myscattering_tot(jfile) = data.scattered.gaiasum + data.scattered.masanasum;
                myscattering_gaia(jfile) = data.scattered.gaiasum;
                myscattering_gaia_psferrpos(jfile) = data.scattered.gaiasum_psferrmax - data.scattered.gaiasum;
                myscattering_gaia_fluxerrpos(jfile) = data.scattered.gaiasum_fluxerrmax - data.scattered.gaiasum;
                myscattering_gaiaerr(jfile) = sqrt(myscattering_gaia_psferrpos(jfile)^2 + myscattering_gaia_fluxerrpos(jfile)^2);
                myscattering_masana(jfile) = data.scattered.masanasum;
                myscattering_masana_psferrpos(jfile) = data.scattered.masanasum_psferrmax - data.scattered.masanasum;
                myscattering_masana_fluxerrpos(jfile) = data.scattered.masanasum_fluxerrmax - data.scattered.masanasum;
                myscattering_masanaerr(jfile) = sqrt(myscattering_masana_psferrpos(jfile)^2 + myscattering_masana_fluxerrpos(jfile)^2);
                myscattering_toterr(jfile) = myscattering_gaiaerr(jfile) + myscattering_masanaerr(jfile);
            end
        end


        % If struct field data.header.bad exists, check if file is good or
        % bad - only mark as good if not bad
        if isfield(data.header,'bad')
            % Only use good images
            if data.header.bad == 0
                isgood(jfile) = 1;
                mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                    mydate(jfile),myfieldnum(jfile),...
                    mytarget{jfile},mysun(jfile),myelong(jfile),...
                    mycrr(jfile),...
                    mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                    mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

                disp(mystring)

                image = data.image.calimage;
                % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
                parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
                parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
            else
                disp(['header.bad WAS BAD on file (',num2str(ifile),'/',num2str(size(wlightfiles,1)),') ',sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name)])
                badcntr = badcntr + 1;
                newest_goodfilesNum = newest_goodfilesNum - 1; %reduce by 1 if there was a bad file
                if( FLG_saveGoodDataFiles == 1 )
                    delete([FLG_saveGoodDataFiles_dir,'/',wlightfiles(ifile).name]) % Remove because it actually was a bad file all along
                end
            end
            % If struct field does not exist, files have not been marked good
            % or bad - assume all good
        else
            isgood(jfile) = 1;
            mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                mydate(jfile),myfieldnum(jfile),...
                mytarget{jfile},mysun(jfile),myelong(jfile),...
                mycrr(jfile),...
                mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

            disp(mystring)

            image = data.image.calimage;
            % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
            parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
            parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
        end
    end
end
disp(['Bad cntr: ',num2str(badcntr)])
disp(['Total files being used for newest data: ',num2str(sum(isgood(1+jfile_offset:numel(wlightfiles)+jfile_offset) & isgoodnewest & ~newest_exclude_logical)),'/',num2str(numel(wlightfiles))])
disp(['isgood (header.bad [also incl. field check and time skip]): ',num2str(sum(isgood(1+jfile_offset:numel(wlightfiles)+jfile_offset))),' | ',...
    'isgoodnewest (good field check): ',num2str(sum(isgoodnewest)),' | ',...
    'newest_exclude (time skip): ',num2str(sum(~newest_exclude_logical))])

%For lauer data files
fprintf('Loading Lauer data \n');
jfile_offset = numoldlightfiles+numnewlightfiles+numnewestlightfiles; %offset to keep ifile in line with other file sets
badcntr = 0;
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
parfor jfile=1+jfile_offset:numel(llightfiles)+jfile_offset
    %If file is for a good field and not excluded, load and save values
    ifile = jfile-jfile_offset; %lets parfor work by letting jfile iterate instead of ifile
    if (isgoodlauer(ifile) == 1) && ~any(lauer_exclude == ifile)

        %Load data files
        dataz = load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(llightfiles,1)));
        data = dataz.data; %get the data out, makes parfor work

        if( FLG_saveGoodDataFiles == 1 )
            [copy_status, copy_msg] = copyfile(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name), [FLG_saveGoodDataFiles_dir,'/',llightfiles(ifile).name]);
            if( copy_status == 0 )
                disp(['WARNING: file ',sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name),' failed to be copied to ',[FLG_saveGoodDataFiles_dir,'/',llightfiles(ifile).name],' failure message follows:'])
                disp(copy_msg)
            end
        end

        %Read in data
        mydate(jfile) = data.header.date_jd-data.header.launch_jd;

        myfieldnum(jfile) = data.header.fieldnum;
        %Assign unique field id based on field number
        unique_field_id(jfile) = find(myfieldnum(jfile) == zones_goodfiles) + zones(zones_indexer)-1; %Counts from 1 to length(goodfiles), like unique_fields
        %         if jfile == 1
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) ~= myfieldnum(jfile-1)
        %             k = k + 1;
        %             unique_field_id(jfile) = unique_fields(k);
        %         elseif myfieldnum(jfile) == myfieldnum(jfile-1)
        %             unique_field_id(jfile) = unique_fields(k);
        %         end
        mydataset(jfile) = zones_indexer; %record the dataset of the data (via the zone, 1 == old, 2 == new, 3 == lauer)

        mytarget{jfile} = data.header.date_cal;

        mysun(jfile) = data.distance.d_sun_NH;

        mydt(jfile) = data.distance.d_NH_target;%./1.496e8;;

        myelong(jfile) = data.coords.sol_elon;

        myb(jfile) = data.coords.galactic(2);

        myexp(jfile) = round(data.header.exptime);

        myref(jfile) = data.ref.mean;

        myeng(jfile) = data.ref.engmean;

        %         myohm(jfile) = data.dgl.ohmmean;

        %         myav(jfile) = 0.4.*data.header.A_V;

        myext(jfile) = data.ext.mean_flux_rat;

        % Reference-corrected mean of masked image
        if want_errmags == 1
            mysig(jfile) = data.err_mags.stats.corrmean;
        else
            mysig(jfile) = data.stats.corrmean;
        end
        mysig_magerr(jfile) = data.err_mags.stats.corrmean;

        mymaskmean(jfile) = data.stats.maskmean;

        maskfrac(jfile) = data.mask.maskfrac;

        % Bias level from header
        %     biaslevl(ifile) = data.astrom.biaslevl;
        % Bias method from header (mean = 1, median = 2) and difference to
        % actual mean/median
        %     if strcmp(data.astrom.biasmthd,'mean of dark column data') == 1
        %         biasmthd(ifile) = 1;
        %         biasdiff(ifile) = data.astrom.biaslevl - mean(data.ref.line);
        %     elseif strcmp(data.astrom.biasmthd,'median of dark column data') == 1
        %         biasmthd(ifile) = 2;
        %         biasdiff(ifile) = data.astrom.biaslevl - median(data.ref.line);
        %     else
        %         fprintf('Absolute panic: bias method not recognized!')
        %     end
        % RSOLAR from header
        rsolar(jfile) = data.astrom.rsolar;

        % Reference-corrected most probable value of masked image
        %     mysig(ifile) = data.stats.corrmostprob;

        % Sum of image brightness of *masked* stars - how bright is the
        % masked-out portion of each image per pixel
        if want_errmags == 1
            starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        else
            starmaskfrac = sum(data.mask.starmask(:))./256.^2;
            mymasked(jfile) = (sum(sum(data.image.calimage.*data.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
        end
        starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
        mymasked_magerr(jfile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);

        mytri(jfile) = data.isl.trimean;
        mytrierr(jfile) = data.isl.trierr;

        if want_errpsf == 1
            mypsfwing(jfile) = data.err_psf.isl.usnowing;
        else
            mypsfwing(jfile) = data.isl.usnowing;
        end
        mypsfwing_psferr(jfile) = data.err_psf.isl.usnowing;

        if sub_tri == 1
            myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);
        elseif sub_tri == 0
            myisl(jfile) = 0;
        end

        pltr_mydgl_planck(jfile) = data.dgl.dglmean_planck;
        pltr_mydglerr_planck(jfile) = data.dgl.dglerr_planck;
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_mydgl_iris(jfile) = data.dgl.dglmean_iris;
        pltr_mydglerr_iris(jfile) = data.dgl.dglerr_iris;
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_mydgl_iris_sfd(jfile) = data.dgl.dglmean_iris_sfd;
        pltr_mydglerr_iris_sfd(jfile) = data.dgl.dglerr_iris_sfd;
        pltr_my100m_planck(jfile) = data.dgl.onehundomean_planck;
        pltr_my100merr_planck(jfile) = data.dgl.ohmstd_planck;
        pltr_my100merr_beta_planck(jfile) = data.dgl.sem_planck_mc_beta;
        pltr_my100merr_temp_planck(jfile) = data.dgl.sem_planck_mc_temp;
        pltr_my100merr_tau_planck(jfile) = data.dgl.sem_planck_mc_tau;
        % Combine temp and beta errors in quad because they're usually equal in magnitude
        pltr_my100merr_combo_planck(jfile) = sqrt(data.dgl.sem_planck_mc_temp^2 + data.dgl.sem_planck_mc_beta^2);
        %     pltr_myohmim_planck(jfile) = data.dgl.ohmim_planck;
        pltr_my100m_iris(jfile) = data.dgl.onehundomean_iris;
        pltr_my100merr_iris(jfile) = std(std(data.dgl.ohmim_iris-0.48));
        %     pltr_myohmim_iris(jfile) = data.dgl.ohmim_iris;
        pltr_my100m_iris_sfd(jfile) = data.dgl.onehundomean_iris_sfd;
        pltr_my100merr_iris_sfd(jfile) = data.dgl.ohmstd_iris_sfd;
        pltr_mydgl_nh(jfile) = data.dgl.ohmmean_nh;
        pltr_mydglerr_nh(jfile) = data.dgl.ohmstd_nh;

        % Calculate mean difference between real cal mean and err sims
        mygalerr(jfile) = mean(data.stats.calmean - data.err_gals.stats.calmean_mc(:));
        mygal_mean(jfile) = mean(data.err_gals.stats.calmean_mc(:));
        mymagerr(jfile) = data.stats.calmean - data.err_mags.stats.calmean;
        mypsferr(jfile) = data.isl.usnowing - data.err_psf.isl.usnowing;

        if want_errmags == 1
            myerr(jfile) = data.err_mags.stats.correrr;
        else
            myerr(jfile) = data.stats.correrr;
        end

        if want_errmags == 1
            mycrr(jfile) = data.err_mags.ref.bias;
        else
            mycrr(jfile) = data.ref.bias;
        end

        if strcmp(flag_method,'new') == 1
            if want_errmags == 1
                myghostdiff(jfile) = data.err_mags.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.err_mags.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.err_mags.ghost.diffusesuberrneg;
            else
                myghostdiff(jfile) = data.ghost.diffusesub;
                myghostdifferrpos(jfile) = data.ghost.diffusesuberrpos;
                myghostdifferrneg(jfile) = data.ghost.diffusesuberrneg;
                %             myghostdiffcomp(jfile) = data.ghost.diffusediff;
                myscattering_tot(jfile) = data.scattered.gaiasum + data.scattered.masanasum;
                myscattering_gaia(jfile) = data.scattered.gaiasum;
                myscattering_gaia_psferrpos(jfile) = data.scattered.gaiasum_psferrmax - data.scattered.gaiasum;
                myscattering_gaia_fluxerrpos(jfile) = data.scattered.gaiasum_fluxerrmax - data.scattered.gaiasum;
                myscattering_gaiaerr(jfile) = sqrt(myscattering_gaia_psferrpos(jfile)^2 + myscattering_gaia_fluxerrpos(jfile)^2);
                myscattering_masana(jfile) = data.scattered.masanasum;
                myscattering_masana_psferrpos(jfile) = data.scattered.masanasum_psferrmax - data.scattered.masanasum;
                myscattering_masana_fluxerrpos(jfile) = data.scattered.masanasum_fluxerrmax - data.scattered.masanasum;
                myscattering_masanaerr(jfile) = sqrt(myscattering_masana_psferrpos(jfile)^2 + myscattering_masana_fluxerrpos(jfile)^2);
                myscattering_toterr(jfile) = myscattering_gaiaerr(jfile) + myscattering_masanaerr(jfile);
            end
        end


        % If struct field data.header.bad exists, check if file is good or
        % bad - only mark as good if not bad
        if isfield(data.header,'bad')
            % Only use good images
            if data.header.bad == 0
                isgood(jfile) = 1;
                mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                    mydate(jfile),myfieldnum(jfile),...
                    mytarget{jfile},mysun(jfile),myelong(jfile),...
                    mycrr(jfile),...
                    mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                    mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

                disp(mystring)

                image = data.image.calimage;
                % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
                parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
                parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
            else
                disp(['header.bad WAS BAD on file (',num2str(ifile),'/',num2str(size(llightfiles,1)),') ',sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name)])
                badcntr = badcntr + 1;
                if( FLG_saveGoodDataFiles == 1 )
                    delete([FLG_saveGoodDataFiles_dir,'/',llightfiles(ifile).name]) % Remove because it actually was a bad file all along
                end
            end
            % If struct field does not exist, files have not been marked good
            % or bad - assume all good
        else
            isgood(jfile) = 1;
            mystring = sprintf('%f, %d, %s, %f, %f, %d, %f, %f, %f, %f',...
                mydate(jfile),myfieldnum(jfile),...
                mytarget{jfile},mysun(jfile),myelong(jfile),...
                mycrr(jfile),...
                mysig(jfile),myisl(jfile),pltr_mydgl_planck(jfile),...
                mysig(jfile)-myisl(jfile)-pltr_mydgl_planck(jfile));

            disp(mystring)

            image = data.image.calimage;
            % save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');
            parsave_image(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),image);

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            % save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
            parsave_nanimage(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),nanimage);
        end
    end
end
disp(['Bad cntr: ',num2str(badcntr)])
disp(['Total files being used for lauer data: ',num2str(sum(isgood(1+jfile_offset:numel(llightfiles)+jfile_offset) & isgoodlauer & ~lauer_exclude_logical)),'/',num2str(numel(llightfiles))])
disp(['isgood (header.bad [also incl. field check and time skip]): ',num2str(sum(isgood(1+jfile_offset:numel(llightfiles)+jfile_offset))),' | ',...
    'isgoodlauer (good field check): ',num2str(sum(isgoodlauer)),' | ',...
    'lauergood_exclude (time skip): ',num2str(sum(~lauer_exclude_logical))])

disp(['Newest|Total good files to be used (including good field (isgoodnewest), time limit (newest_exclude), isgood (header.bad check)): ',num2str(newest_goodfilesNum)])
% Select only data marked as good

pltr_thissig = mysig(isgood)-myisl(isgood)-myghostdiff(isgood)-myscattering_tot(isgood);

isgood_masked = mymasked(isgood);
isgood_corrmean = mysig(isgood);
isgood_maskmean = mymaskmean(isgood);
isgood_tri = mytri(isgood);
isgood_psfwing = mypsfwing(isgood);
isgood_isl = myisl(isgood);
isgood_dgl = pltr_mydgl_planck(isgood);
isgood_trierr = mytrierr(isgood);
isgood_dglerr = pltr_mydglerr_planck(isgood);
if strcmp(flag_method,'new') == 1
    isgood_diffghost = myghostdiff(isgood);
    isgood_diffghosterrpos = myghostdifferrpos(isgood);
    isgood_diffghosterrneg = myghostdifferrneg(isgood);
    if want_errmags == 0 && want_errpsf == 0
        isgood_scattering_tot = myscattering_tot(isgood);
        isgood_scattering_gaia = myscattering_gaia(isgood);
        isgood_scattering_gaia_psferrpos = myscattering_gaia_psferrpos(isgood);
        isgood_scattering_gaia_fluxerrpos = myscattering_gaia_fluxerrpos(isgood);
        isgood_scattering_masana = myscattering_masana(isgood);
        isgood_scattering_masana_psferrpos = myscattering_masana_psferrpos(isgood);
        isgood_scattering_masana_fluxerrpos = myscattering_masana_fluxerrpos(isgood);
    end
end
isgood_ext = myext(isgood);
isgood_field = unique_field_id(isgood);
isgood_b = myb(isgood);

% Calculate per-field values
goodfiles_ditch = false(length(goodfiles),1);
for i = 1:length(goodfiles)
    k_field = isgood_field == unique_fields(i);
    if( sum(k_field) > 1 )
        field_masked = mean(isgood_masked(k_field));
        field_maskmean = mean(isgood_maskmean(k_field));
        field_tri = mean(isgood_tri(k_field));
        field_psfwing = mean(isgood_psfwing(k_field));
        field_corrmean = mean(isgood_corrmean(k_field));
        field_isl = mean(isgood_isl(k_field));
        field_dgl = mean(isgood_dgl(k_field));
        field_trierr = mean(isgood_trierr(k_field));
        field_dglerr = mean(isgood_dglerr(k_field));
        field_ext = mean(isgood_ext(k_field));
        if strcmp(flag_method,'new') == 1
            field_diffghost = mean(isgood_diffghost(k_field));
            field_diffghosterrpos = mean(isgood_diffghosterrpos(k_field));
            field_diffghosterrneg = mean(isgood_diffghosterrneg(k_field));
            if want_errmags == 0 && want_errpsf == 0
                field_scattering_tot = mean(isgood_scattering_tot(k_field));
                field_scattering_gaia = mean(isgood_scattering_gaia(k_field));
                field_scattering_gaia_psferrpos = mean(isgood_scattering_gaia_psferrpos(k_field));
                field_scattering_gaia_fluxerrpos = mean(isgood_scattering_gaia_fluxerrpos(k_field));
                field_scattering_masana = mean(isgood_scattering_masana(k_field));
                field_scattering_masana_psferrpos = mean(isgood_scattering_masana_psferrpos(k_field));
                field_scattering_masana_fluxerrpos = mean(isgood_scattering_masana_fluxerrpos(k_field));
            else
                field_scattering_tot = 0;
                field_scattering_gaia = 0;
                field_scattering_gaia_psferrpos = 0;
                field_scattering_gaia_fluxerrpos = 0;
                field_scattering_masana = 0;
                field_scattering_masana_psferrpos = 0;
                field_scattering_masana_fluxerrpos = 0;
            end
        else
            field_diffghost = 0;
            field_diffghosterrpos = 0;
            field_diffghosterrneg = 0;
            field_scattering = 0;
        end
        field_b = mean(isgood_b(k_field));
        disp(['Field #',num2str(goodfiles(i)),': b = ' ,num2str(field_b),...
            ' | masked mean = ',num2str(field_masked),' | corr mean = ',num2str(field_corrmean),...
            ' | tri mean = ',num2str(field_tri),' | tri err = ',num2str(field_trierr),...
            ' |  psfwing mean = ',num2str(field_psfwing),' | ISL mean = ',num2str(field_isl),...
            ' | Planck DGL mean = ',num2str(field_dgl),' | Planck DGL err = ',num2str(field_dglerr),...
            ' | Scattering mean = ',num2str(field_scattering_tot),' | Scattering Gaia = ',num2str(field_scattering_gaia),' | Scattering Masana = ',num2str(field_scattering_masana),...
            ' | Scattering Masana PSF err = ',num2str(field_scattering_masana_psferrpos),' | Scattering Masana Flux err = ',num2str(field_scattering_masana_fluxerrpos),...
            ' | Scattering Gaia PSF err = ',num2str(field_scattering_gaia_psferrpos),' | Scattering Gaia Flux err = ',num2str(field_scattering_gaia_fluxerrpos),...
            ' | Diff Ghost mean = ',num2str(field_diffghost),' | Diff Ghost err pos = ',num2str(field_diffghosterrpos),' | Diff Ghost err neg = ',num2str(field_diffghosterrneg),...
            ' | Ext = ', num2str(field_ext)])
    else
        %if k_field is empty then no goodfiles in the field
        %remove the field from the listing, can't do it mid-loop though
        goodfiles_ditch(i) = 1; %set to ditch
    end
end
%now remove empty fields
goodfiles(goodfiles_ditch) = [];
unique_fields(goodfiles_ditch) = [];

figcnt = 1;

figure(figcnt); clf
plot(mysun(isgood),mysig(isgood)-myisl(isgood)-pltr_mydgl_planck(isgood),'o')
xlabel('Solar Distance')
ylabel('EBL')
figcnt = figcnt + 1;

if strcmp(flag_method,'new') == 1
    mysun_isgood = mysun(isgood);
    myghostdiff_isgood = myghostdiff(isgood);
    
    figure(figcnt); clf
    hold on;
    pHolder = gobjects(length(unique_fields),1);
    pNamer = cell(length(unique_fields),1);
    for jk = 1:length(unique_fields)
        field_curr = find(isgood_field == unique_fields(jk));
        pHolder(jk) = plot(mysun_isgood(field_curr),myghostdiff_isgood(field_curr),'o','MarkerFaceColor',PLOT_color{jk},'color',PLOT_color{jk});
        pNamer{jk} = ['Field ',num2str(jk)];
    end
    legend(pHolder,pNamer,'NumColumns',2) %'Location','northwest'
%     hold on;
%     errorbar(mysun(isgood),myghostdiff(isgood), myghostdifferrneg(isgood), myghostdifferrpos(isgood), 'LineStyle','none','Color','k');
    xlabel('Solar Distance')
    ylabel('Summed Ghost Intensity [nW m^{-2} sr^{-1}]')
    figcnt = figcnt + 1;

    %     figure(figcnt); clf
    %     plot(mysun(isgood),myghostdiffcomp(isgood),'o')
    %     yline(mean(myghostdiffcomp(isgood)))
    %     xlabel('Solar Distance')
    %     ylabel('(Predicted Summed Ghost Intensity - ISL Estimation) [nW m^{-2} sr^{-1}]')
%     figcnt = figcnt + 1;

end

% Preallocate per-field values
mydist = zeros(numel(goodfiles),1);
mygal = zeros(numel(goodfiles),1);
mysubmen = zeros(numel(goodfiles),1);
mysuberr = zeros(numel(goodfiles),1);
mymean = zeros(numel(goodfiles),1);
myunc = zeros(numel(goodfiles),1);
mysem = zeros(numel(goodfiles),1);
% myavp = zeros(numel(goodfiles),1);
myohmp = zeros(numel(goodfiles),1);
myb_goodfiles = zeros(numel(goodfiles),1);
mysig_goodfiles = zeros(numel(goodfiles),1);
mysig_magerr_goodfiles = zeros(numel(goodfiles),1);
mymaskmean_goodfiles = zeros(numel(goodfiles),1);
mymasked_goodfiles = zeros(numel(goodfiles),1);
mymasked_magerr_goodfiles = zeros(numel(goodfiles),1);
mytri_goodfiles = zeros(numel(goodfiles),1);
mypsfwing_goodfiles = zeros(numel(goodfiles),1);
mypsfwing_psferr_goodfiles = zeros(numel(goodfiles),1);
myisl_goodfiles = zeros(numel(goodfiles),1);
myghostdiff_goodfiles = zeros(numel(goodfiles),1);
myscattering_tot_goodfiles = zeros(numel(goodfiles),1);
myscattering_masana_goodfiles = zeros(numel(goodfiles),1);
myscattering_gaia_goodfiles = zeros(numel(goodfiles),1);
mygal_mean_goodfiles = zeros(numel(goodfiles),1);
pltr_thissig_goodfiles = zeros(numel(goodfiles),1);
pltr_mydgl_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_mydglerr_planck_goodfiles = zeros(numel(goodfiles),1);
% pltr_myohmim_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_mydgl_iris_goodfiles = zeros(numel(goodfiles),1);
pltr_mydglerr_iris_goodfiles = zeros(numel(goodfiles),1);
% pltr_myohmim_iris_goodfiles = zeros(numel(goodfiles),1);
pltr_mydgl_iris_sfd_goodfiles = zeros(numel(goodfiles),1);
pltr_mydglerr_iris_sfd_goodfiles = zeros(numel(goodfiles),1);
pltr_my100m_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_beta_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_temp_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_tau_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_combo_planck_goodfiles = zeros(numel(goodfiles),1);
pltr_my100m_iris_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_iris_goodfiles = zeros(numel(goodfiles),1);
pltr_my100m_iris_sfd_goodfiles = zeros(numel(goodfiles),1);
pltr_my100merr_iris_sfd_goodfiles = zeros(numel(goodfiles),1);
pltr_mydgl_nh_goodfiles = zeros(numel(goodfiles),1);
pltr_mydglerr_nh_goodfiles = zeros(numel(goodfiles),1);
pltr_sem = zeros(numel(goodfiles),1);
pltr_unc = zeros(numel(goodfiles),1);
pltr_mydataset = zeros(numel(goodfiles),1);
pltr_mydatasetTriplets = zeros(numel(goodfiles),3);
mygalerr_goodfiles = zeros(numel(goodfiles),1);
myext_goodfiles = zeros(numel(goodfiles),1);
mymagerr_goodfiles = zeros(numel(goodfiles),1);
mypsferr_goodfiles = zeros(numel(goodfiles),1);
mytrierr_goodfiles = zeros(numel(goodfiles),1);
myscattering_toterr_goodfiles = zeros(numel(goodfiles),1);
myscattering_masanaerr_goodfiles = zeros(numel(goodfiles),1);
myscattering_gaiaerr_goodfiles = zeros(numel(goodfiles),1);
myghostdifferrpos_goodfiles = zeros(numel(goodfiles),1);
myghostdifferrneg_goodfiles = zeros(numel(goodfiles),1);
mytoterr_pos = zeros(numel(goodfiles),1);
mytoterr_neg = zeros(numel(goodfiles),1);
myfieldnum_goodfiles = zeros(numel(goodfiles),1);
myfieldexp_goodfiles = zeros(numel(goodfiles),1);
maskfrac_goodfiles = zeros(numel(goodfiles),1);
solelon_goodfiles = zeros(numel(goodfiles),1);
mysun_goodfiles = zeros(numel(goodfiles),1);

% Calculate per-field values
for ifield=1:numel(goodfiles)

    whpl = unique_field_id == unique_fields(ifield);
    mydist(ifield) = mean(mysun(whpl & isgood));
    mygal(ifield) = mean(myb(whpl & isgood));
    myfieldnum_goodfiles(ifield) = mean(myfieldnum(whpl & isgood));
    myfieldexp_goodfiles(ifield) = length(myb(whpl & isgood));
    if strcmp(flag_method,'new') == 1
        thissig_planck = mysig(whpl & isgood)-myisl(whpl & isgood)-pltr_mydgl_planck(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        thissig_iris = mysig(whpl & isgood)-myisl(whpl & isgood)-pltr_mydgl_iris(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        thissig_iris_sfd = mysig(whpl & isgood)-myisl(whpl & isgood)-pltr_mydgl_iris_sfd(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        pltr_totsig_goodfiles = mysig(whpl & isgood)-myisl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        pltr_thissig_goodfiles(ifield) = mean(mysig(whpl & isgood)-myisl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood));
    else
        thissig = mysig(whpl & isgood)-myisl(whpl & isgood)-pltr_mydgl_planck(whpl & isgood);
    end
    mysig_goodfiles(ifield) = mean(mysig(whpl & isgood));
    mysig_magerr_goodfiles(ifield) = mean(mysig_magerr(whpl & isgood));
    mymaskmean_goodfiles(ifield) = mean(mymaskmean(whpl & isgood));
    mymasked_goodfiles(ifield) = mean(mymasked(whpl & isgood));
    mymasked_magerr_goodfiles(ifield) = mean(mymasked_magerr(whpl & isgood));
    mytri_goodfiles(ifield) = mean(mytri(whpl & isgood));
    mypsfwing_goodfiles(ifield) = mean(mypsfwing(whpl & isgood));
    mypsfwing_psferr_goodfiles(ifield) = mean(mypsfwing_psferr(whpl & isgood));
    myisl_goodfiles(ifield) = mean(myisl(whpl & isgood));
    myghostdiff_goodfiles(ifield) = mean(myghostdiff(whpl & isgood));
    myscattering_tot_goodfiles(ifield) = mean(myscattering_tot(whpl & isgood));
    myscattering_masana_goodfiles(ifield) = mean(myscattering_masana(whpl & isgood));
    myscattering_gaia_goodfiles(ifield) = mean(myscattering_gaia(whpl & isgood));
    mygal_mean_goodfiles(ifield) = mean(mygal_mean(whpl & isgood));
    maskfrac_goodfiles(ifield) = mean(maskfrac(whpl & isgood));
    solelon_goodfiles(ifield) = mean(myelong(whpl & isgood));
    mysun_goodfiles(ifield) = mean(mysun(whpl & isgood));

    pltr_mydgl_planck_goodfiles(ifield) = mean(pltr_mydgl_planck(whpl & isgood));
    pltr_mydglerr_planck_goodfiles(ifield) = mean(pltr_mydglerr_planck(whpl & isgood));
    %     pltr_myohmim_planck_goodfiles(ifield) = pltr_myohmim_planck(whpl);
    pltr_mydgl_iris_goodfiles(ifield) = mean(pltr_mydgl_iris(whpl & isgood));
    pltr_mydglerr_iris_goodfiles(ifield) = mean(pltr_mydglerr_iris(whpl & isgood));
    %     pltr_myohmim_iris_goodfiles(ifield) = pltr_myohmim_iris(whpl);
    pltr_mydgl_iris_sfd_goodfiles(ifield) = mean(pltr_mydgl_iris_sfd(whpl & isgood));
    pltr_mydglerr_iris_sfd_goodfiles(ifield) = mean(pltr_mydglerr_iris_sfd(whpl & isgood));
    pltr_my100m_planck_goodfiles(ifield) = mean(pltr_my100m_planck(whpl & isgood));
    pltr_my100merr_planck_goodfiles(ifield) = mean(pltr_my100merr_planck(whpl & isgood));
    pltr_my100merr_beta_planck_goodfiles(ifield) = mean(pltr_my100merr_beta_planck(whpl & isgood));
    pltr_my100merr_temp_planck_goodfiles(ifield) = mean(pltr_my100merr_temp_planck(whpl & isgood));
    pltr_my100merr_tau_planck_goodfiles(ifield) = mean(pltr_my100merr_tau_planck(whpl & isgood));
    pltr_my100merr_combo_planck_goodfiles(ifield) = mean(pltr_my100merr_combo_planck(whpl & isgood));
    pltr_my100m_iris_goodfiles(ifield) = mean(pltr_my100m_iris(whpl & isgood));
    pltr_my100merr_iris_goodfiles(ifield) = mean(pltr_my100merr_iris(whpl & isgood));
    pltr_my100m_iris_sfd_goodfiles(ifield) = mean(pltr_my100m_iris_sfd(whpl & isgood));
    pltr_my100merr_iris_sfd_goodfiles(ifield) = mean(pltr_my100merr_iris_sfd(whpl & isgood));
    pltr_mydgl_nh_goodfiles(ifield) = mean(pltr_mydgl_nh(whpl & isgood));
    pltr_mydglerr_nh_goodfiles(ifield) = mean(pltr_mydglerr_nh(whpl & isgood));

    myb_goodfiles(ifield) = mean(myb(whpl & isgood));
    mygalerr_goodfiles(ifield) = mean(mygalerr(whpl&isgood));
    myext_goodfiles(ifield) = mean(myext(whpl&isgood));
    mymagerr_goodfiles(ifield) = mean(mymagerr(whpl&isgood));
    mypsferr_goodfiles(ifield) = mean(mypsferr(whpl&isgood));
    mytrierr_goodfiles(ifield) = mean(mytrierr(whpl&isgood));
    myscattering_toterr_goodfiles(ifield) = mean(myscattering_toterr(whpl&isgood));
    myscattering_masanaerr_goodfiles(ifield) = mean(myscattering_masanaerr(whpl&isgood));
    myscattering_gaiaerr_goodfiles(ifield) = mean(myscattering_gaiaerr(whpl&isgood));
    myghostdifferrpos_goodfiles(ifield) = mean(myghostdifferrpos(whpl&isgood));
    myghostdifferrneg_goodfiles(ifield) = mean(myghostdifferrneg(whpl&isgood));

    thiserr = myerr(whpl & isgood);
    mysubmen(ifield) = sum(mysig(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    mysuberr(ifield) = std(mysig(whpl & isgood));
    mymean_planck(ifield) = sum(thissig_planck./thiserr.^2)./sum(1./thiserr.^2);
    myunc_planck(ifield) = std(thissig_planck);
    mysem_planck(ifield) = std(thissig_planck)/sqrt(length(thissig_planck));
    mymean_iris(ifield) = sum(thissig_iris./thiserr.^2)./sum(1./thiserr.^2);
    myunc_iris(ifield) = std(thissig_iris);
    mysem_iris(ifield) = std(thissig_iris)/sqrt(length(thissig_iris));
    mymean_iris_sfd(ifield) = sum(thissig_iris_sfd./thiserr.^2)./sum(1./thiserr.^2);
    myunc_iris_sfd(ifield) = std(thissig_iris_sfd);
    mysem_iris_sfd(ifield) = std(thissig_iris_sfd)/sqrt(length(thissig_iris_sfd));
    pltr_sem(ifield) = std(pltr_totsig_goodfiles)/sqrt(length(pltr_totsig_goodfiles));
    pltr_unc(ifield) = std(pltr_totsig_goodfiles);
    %     myavp(ifield) = sum(myav(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    %     myohmp(ifield) = sum(myohm(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    pltr_mydataset(ifield) = median(mydataset(whpl)); %for plotting per dataset color
%     mytoterr_pos(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
%         mytrierr_goodfiles(ifield)^2 + myghostdifferrpos_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + pltr_sem(ifield)^2);
%     mytoterr_neg(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
%         mytrierr_goodfiles(ifield)^2 + myghostdifferrneg_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + pltr_sem(ifield)^2);
    
    % Use std (myunc) instead of sem (mysem) for individual statistical
    % errors
    if sub_tri == 1
        mytoterr_pos(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
        mytrierr_goodfiles(ifield)^2 + myghostdifferrpos_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + pltr_unc(ifield)^2);
    elseif sub_tri == 0
        mytoterr_pos(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + ...
        myghostdifferrpos_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + pltr_unc(ifield)^2);
    end
    mytoterr_neg(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
        mytrierr_goodfiles(ifield)^2 + myghostdifferrneg_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + pltr_unc(ifield)^2);
end
pltr_mydatasetTriplets(pltr_mydataset == 1,:) = repmat([0, 0.4470, 0.7410],sum(pltr_mydataset == 1),1); %set consistent color triplets
pltr_mydatasetTriplets(pltr_mydataset == 2,:) = repmat([0.8500, 0.3250, 0.0980],sum(pltr_mydataset == 2),1);
pltr_mydatasetTriplets(pltr_mydataset == 3,:) = repmat([0.9290, 0.6940, 0.1250],sum(pltr_mydataset == 3),1);

% figure(figcnt); clf
% plot(mydist,mysubmen,'o')
% hold on
% errorbar(mydist,mysubmen,mysuberr,'o');
% xlabel('Solar Distance')
% ylabel('Error-Weighted Image Mean')
% figcnt = figcnt + 1;

% mysubmen

% Direct COB sub - Planck
figure(figcnt); clf
FLG_perDatasetColoring = false; %false plots default (same color), true plots colors per dataset (old/new/lar)
%plot(mydist,mymean,'.','MarkerSize',0);
%hold on
if ~FLG_perDatasetColoring
    errorbar(mydist,mymean_planck,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
else
    hold on;
    for j = 1:length(mydist)
        errorbar(mydist(j),mymean_planck(j),myunc_planck(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
end
hold on;

supermean = sum(mymean_planck./myunc_planck.^2)./sum(1./myunc_planck.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./myunc_planck.^2)) %old way was myunc instead of mysem
supersub = sum(mysubmen./myunc_planck.^2)./sum(1./myunc_planck.^2)

plot(xlim,[supermean,supermean],'Color',[0, 0.4470, 0.7410],'LineWidth',2)
plot(xlim,[supermean+superunc,supermean+superunc],'Color',[0, 0.4470, 0.7410],'LineStyle','-.','LineWidth',1)
plot(xlim,[supermean-superunc,supermean-superunc],'Color',[0, 0.4470, 0.7410],'LineStyle','-.','LineWidth',1)
plot(xlim,[0,0],'k--','LineWidth',2)
ylim([-30,25]);
xlim([5,50]);
xlabel('Heliocentric Distance [AU]')
ylabel('Direct-Subtraction COB [nW m^{-2} sr^{-1}]')
title('Planck')
figcnt = figcnt + 1;

% Direct COB sub - IRIS
figure(figcnt); clf
FLG_perDatasetColoring = false; %false plots default (same color), true plots colors per dataset (old/new/lar)
%plot(mydist,mymean,'.','MarkerSize',0);
%hold on
if ~FLG_perDatasetColoring
    errorbar(mydist,mymean_iris,myunc_iris,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
else
    hold on;
    for j = 1:length(mydist)
        errorbar(mydist(j),mymean_iris(j),myunc_iris(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
end
hold on;

supermean = sum(mymean_iris./myunc_iris.^2)./sum(1./myunc_iris.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./myunc_iris.^2)) %old way was myunc instead of mysem
supersub = sum(mysubmen./myunc_iris.^2)./sum(1./myunc_iris.^2)

plot(xlim,[supermean,supermean],'Color',[0, 0.4470, 0.7410],'LineWidth',2)
plot(xlim,[supermean+superunc,supermean+superunc],'Color',[0, 0.4470, 0.7410],'LineStyle','-.','LineWidth',1)
plot(xlim,[supermean-superunc,supermean-superunc],'Color',[0, 0.4470, 0.7410],'LineStyle','-.','LineWidth',1)
plot(xlim,[0,0],'k--','LineWidth',2)
ylim([-30,25]);
xlim([5,50]);
xlabel('Heliocentric Distance [AU]')
ylabel('Direct-Subtraction COB [nW m^{-2} sr^{-1}]')
title('IRIS')
figcnt = figcnt + 1;

% Direct COB sub - IRIS/SFD
figure(figcnt); clf
FLG_perDatasetColoring = false; %false plots default (same color), true plots colors per dataset (old/new/lar)
%plot(mydist,mymean,'.','MarkerSize',0);
%hold on
if ~FLG_perDatasetColoring
    errorbar(mydist,mymean_iris_sfd,myunc_iris_sfd,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
else
    hold on;
    for j = 1:length(mydist)
        errorbar(mydist(j),mymean_iris_sfd(j),myunc_iris_sfd(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
end
hold on;

supermean = sum(mymean_iris_sfd./myunc_iris_sfd.^2)./sum(1./myunc_iris_sfd.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./myunc_iris_sfd.^2)) %old way was myunc instead of mysem
supersub = sum(mysubmen./myunc_iris_sfd.^2)./sum(1./myunc_iris_sfd.^2)

plot(xlim,[supermean,supermean],'Color',[0, 0.4470, 0.7410],'LineWidth',2)
plot(xlim,[supermean+superunc,supermean+superunc],'Color',[0, 0.4470, 0.7410],'LineStyle','-.','LineWidth',1)
plot(xlim,[supermean-superunc,supermean-superunc],'Color',[0, 0.4470, 0.7410],'LineStyle','-.','LineWidth',1)
plot(xlim,[0,0],'k--','LineWidth',2)
ylim([-30,25]);
xlim([5,50]);
xlabel('Heliocentric Distance [AU]')
ylabel('Direct-Subtraction COB [nW m^{-2} sr^{-1}]')
title('IRIS/SFD')
figcnt = figcnt + 1;

%=========================combo of Planck, IRIS, IRIS/SFD=========================
figure(figcnt); clf
hold on;
FLG_perDatasetColoring = false; %false plots default (same color), true plots colors per dataset (old/new/lar)

pHolder = gobjects(3,1);
pNamer = {'Planck','IRIS','IRIS/SFD'};
pOffset = [0,2,4];
xlim_r = [5,50];

if ~FLG_perDatasetColoring
    %Planck
    pHolder(1) = errorbar(mydist+pOffset(1),mymean_planck,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
    %IRIS
    pHolder(2) = errorbar(mydist+pOffset(2),mymean_iris,myunc_iris,'.','MarkerSize',20,'MarkerEdge',[0.4940, 0.1840, 0.5560],'LineStyle','none','Color',[0.4940, 0.1840, 0.5560]); %old error bars were myunc, based on std not sem
    %IRIS/SFD
    pHolder(3) = errorbar(mydist+pOffset(3),mymean_iris_sfd,myunc_iris_sfd,'.','MarkerSize',20,'MarkerEdge',[0.4660, 0.6740, 0.1880],'LineStyle','none','Color',[0.4660, 0.6740, 0.1880]); %old error bars were myunc, based on std not sem
else
    %Planck
    for j = 1:length(mydist)
        pHolder(1) = errorbar(mydist(j)+pOffset(1),mymean_planck(j),myunc_planck(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
    %IRIS
    for j = 1:length(mydist)
        pHolder(2) = errorbar(mydist(j)+pOffset(2),mymean_iris(j),myunc_iris(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
    %IRIS/SFD
    for j = 1:length(mydist)
        pHolder(3) = errorbar(mydist(j)+pOffset(3),mymean_iris_sfd(j),myunc_iris_sfd(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
end

%print some stuff?
disp('Planck')
supermean = sum(mymean_planck./myunc_planck.^2)./sum(1./myunc_planck.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./myunc_planck.^2)) %old way was myunc instead of mysem
supersub = sum(mysubmen./myunc_planck.^2)./sum(1./myunc_planck.^2)
fill([xlim_r,flip(xlim_r)],[supermean+superunc,supermean+superunc,supermean-superunc,supermean-superunc],[0.8500, 0.3250, 0.0980],'LineStyle','none');
alpha(.5);
disp('IRIS')
supermean = sum(mymean_iris./myunc_iris.^2)./sum(1./myunc_iris.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./myunc_iris.^2)) %old way was myunc instead of mysem
supersub = sum(mysubmen./myunc_iris.^2)./sum(1./myunc_iris.^2)
fill([xlim_r,flip(xlim_r)],[supermean+superunc,supermean+superunc,supermean-superunc,supermean-superunc],[0.4940, 0.1840, 0.5560],'LineStyle','none');
alpha(.5);
disp('IRIS/SFD')
supermean = sum(mymean_iris_sfd./myunc_iris_sfd.^2)./sum(1./myunc_iris_sfd.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./myunc_iris_sfd.^2)) %old way was myunc instead of mysem
supersub = sum(mysubmen./myunc_iris_sfd.^2)./sum(1./myunc_iris_sfd.^2)
fill([xlim_r,flip(xlim_r)],[supermean+superunc,supermean+superunc,supermean-superunc,supermean-superunc],[0.4660, 0.6740, 0.1880],'LineStyle','none');
alpha(.5);

plot(xlim,[0,0],'k--','LineWidth',2)
ylim([-30,25]);
xlim(xlim_r);
xlabel('Heliocentric Distance [AU]')
ylabel('Direct-Subtraction COB [nW m^{-2} sr^{-1}]')
legend(pHolder,pNamer);
figcnt = figcnt + 1;


figure(figcnt); clf
me = errorbar(abs(mygal),mymean_planck,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
hold on;

z_pts = [15.4,18.1,-6.3,29.2];
z_err = [14.4,26.2,9.1,20.5];
z_err_sem = [14.4/sqrt(10),26.2/sqrt(10),9.1/sqrt(3),20.5/sqrt(3)];
z_b = [85.74,28.41,57.69,62.03];

zem = errorbar(abs(z_b),z_pts,z_err,'.','MarkerSize',20,'MarkerEdge',[0, 0.4470, 0.7410],'LineStyle','none','Color',[0, 0.4470, 0.7410]); %old error bars were myunc, based on std not sem
% supermean = sum(mymean./myunc.^2)./sum(1./myunc.^2) %old way was myunc instead of mysem
% superunc = 1./sqrt(sum(1./myunc.^2)) %old way was myunc instead of mysem
%
% superav = sum(myavp./myunc.^2)./sum(1./myunc.^2)
%
% supersub = sum(mysubmen./myunc.^2)./sum(1./myunc.^2)
%
% plot(xlim,[supermean,supermean],'Color',[0, 0.4470, 0.7410])
% plot(xlim,[supermean+superunc,supermean+superunc],'Color',[0, 0.4470, 0.7410],'LineStyle',':')
% plot(xlim,[supermean-superunc,supermean-superunc],'Color',[0, 0.4470, 0.7410],'LineStyle',':')
% plot(xlim,[0,0],'k:')
legend([me,zem],{'Symons Fields','Zemcov Fields'},'Location','southwest');
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('COB [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;

% Plot DGL vs. EBL + DGL for all images
% figure(figcnt); clf
% p1 = scatter(pltr_mydgl_planck(isgood),pltr_thissig,'o');
% hold on
% % p1err = errorbar(pltr_thissig,pltr_mydgl_planck,pltr_mydglerr_planck,'o');
% i1 = scatter(pltr_mydgl_iris(isgood),pltr_thissig,'g');
% % i1err = errorbar(pltr_thissig,pltr_mydgl_iris,pltr_mydglerr_iris,'m');
% ylabel('EBL + DGL (nW/m^2/sr)')
% xlabel('DGL (nW/m^2/sr)')
% legend([p1,i1],{'Planck','Iris'},'Location','northwest');
% figcnt = figcnt + 1;

% Plot DGL vs. EBL + DGL per field
% figure(figcnt); clf
% p1 = scatter(pltr_mydgl_planck_goodfiles,pltr_thissig_goodfiles,'o');
% hold on
% % p1err = errorbar(pltr_thissig,pltr_mydgl_planck,pltr_mydglerr_planck,'o');
% i1 = scatter(pltr_mydgl_iris_goodfiles,pltr_thissig_goodfiles,'g');
% % i1err = errorbar(pltr_thissig,pltr_mydgl_iris,pltr_mydglerr_iris,'m');
% ylabel('EBL + DGL (nW/m^2/sr)')
% xlabel('DGL (nW/m^2/sr)')
% legend([p1,i1],{'Planck','Iris'},'Location','northwest');
% figcnt = figcnt + 1;

% Plot 100 micron emission vs. EBL + DGL with fit that includes x and y errors
figure(figcnt); clf
FLG_perDatasetColoring = true; %false plots colors per Planck/Iris/Iris+SDF, true plots colors per dataset (old/new/lar)

% IRIS error
iris_err = (0.06)/sqrt((1.13*(4.3)^2)/((17.4^2))); % rms noise (0.06 MJy/sr) converted from per iris beam to per lorri image
iris_err_fields = ones(length(goodfiles),1)*iris_err;
% Scaling factor for 100 micron emission with galactic latitude (free
% parameter A not included b/c part of b(lambda) - OR IS IT??)
dl = dglparams.A.*(1 - 1.1 .* (dglparams.g(1)).*sqrt(sin(abs(myb_goodfiles).*pi./180)));
% Calculate propagated error on d(b) due to error on g
dberr = (dglparams.g(2).^2 .* (dglparams.A.*1.1 .*sqrt(sin(abs(myb_goodfiles).*pi./180))).^2); %Have added A into this
% Calculate propagated error on 100m*d(b) due to error on 100m and error on d(b)
planck_prop_err = sqrt(((pltr_my100merr_combo_planck_goodfiles.*dl).^2) + (dberr.*(pltr_my100m_planck_goodfiles.^2)));
iris_prop_err = sqrt(((iris_err_fields.*dl).^2) + (dberr.*(pltr_my100m_iris_goodfiles.^2)));
iris_sfd_prop_err = sqrt(((iris_err_fields.*dl).^2) + (dberr.*(pltr_my100m_iris_sfd_goodfiles.^2)));

% make fit
% planck
[b_planck, m_planck, sigm_planck, sigb_planck, chi2, q] = fitexy(pltr_my100m_planck_goodfiles.*dl, pltr_thissig_goodfiles, planck_prop_err, mytoterr_pos);
% iris
[b_iris, m_iris, sigm_iris, sigb_iris, chi2, q] = fitexy(pltr_my100m_iris_goodfiles.*dl, pltr_thissig_goodfiles, iris_prop_err, mytoterr_pos);
% iris/sfd
[b_iris_sfd, m_iris_sfd, sigm_iris_sfd, sigb_iris_sfd, chi2, q] = fitexy(pltr_my100m_iris_sfd_goodfiles.*dl, pltr_thissig_goodfiles, iris_sfd_prop_err, mytoterr_pos);

if ~FLG_perDatasetColoring
    p1 = scatter(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,'b');
    hold on
    i1 = scatter(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,'g');
    i2 = scatter(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,'m');
else
    p1 = scatter(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,200,pltr_mydatasetTriplets,'filled');
    hold on
    i1 = scatter(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,200,pltr_mydatasetTriplets,'filled');
    i2 = scatter(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,200,pltr_mydatasetTriplets,'filled');
end
p1errx = errorbar(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,pltr_my100merr_combo_planck_goodfiles.*dl,'horizontal', 'LineStyle','none','Color', 'b');
p1erry = errorbar(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,mytoterr_pos,'vertical', 'LineStyle','none','Color', 'b');

i1errx = errorbar(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,iris_err_fields.*dl,'horizontal', 'LineStyle','none','Color', 'g');
i1erry = errorbar(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,mytoterr_pos,'vertical', 'LineStyle','none','Color', 'g');

i2errx = errorbar(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,iris_err_fields.*dl,'horizontal', 'LineStyle','none','Color', 'm');
i2erry = errorbar(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,mytoterr_pos,'vertical', 'LineStyle','none','Color', 'm');

% apply fit
% fit_x_planck = linspace(min(pltr_my100m_planck_goodfiles),max(pltr_my100m_planck_goodfiles));
fit_x_planck = linspace(0,max(pltr_my100m_planck_goodfiles.*dl));
fit_y_planck = (m_planck*fit_x_planck + b_planck);
% fit_x_iris = linspace(min(pltr_my100m_iris_goodfiles),max(pltr_my100m_iris_goodfiles));
fit_x_iris = linspace(0,max(pltr_my100m_iris_goodfiles.*dl));
fit_y_iris = (m_iris*fit_x_iris + b_iris);
% fit_x_iris_sfd = linspace(min(pltr_my100m_iris_sfd_goodfiles),max(pltr_my100m_iris_sfd_goodfiles));
fit_x_iris_sfd = linspace(0,max(pltr_my100m_iris_sfd_goodfiles.*dl));
fit_y_iris_sfd = (m_iris_sfd*fit_x_iris_sfd + b_iris_sfd);

% plot fit
fit_planck = plot(fit_x_planck,fit_y_planck,'b');
fit_iris = plot(fit_x_iris,fit_y_iris,'g');
fit_iris_sfd = plot(fit_x_iris_sfd,fit_y_iris_sfd,'m');

ylabel('EBL + DGL (nW/m^2/sr)')
xlabel('100 \mum (MJy/sr) * (1-1.1*g*sin(b)^{1/2})')
legend([p1,i1,i2],{'Planck','Iris','Iris + SFD'},'Location','northwest');
figcnt = figcnt + 1;

% Plot NHI vs. EBL + DGL
% figure(figcnt); clf
% %
% % % HI4PI error
% hi4pi_err = ((2.3e18)/5)/sqrt((1.13*(16.2)^2)/((17.4^2))); % 5-sig rms sensitivity (2.3e18 cm^-2) converted to 1-sig and from per hi4pi beam to per lorri image
% hi4pi_err_fields = ones(length(goodfiles),1)*hi4pi_err;
% n1 = scatter(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,'k');
% hold on
% n1errx = errorbar(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,hi4pi_err_fields,'horizontal', 'LineStyle','none','Color', 'k');
% n1erry = errorbar(pltr_mydgl_nh_goodfiles,pltr_thisiris_err_fields.*dlsig_goodfiles,mytoterr_pos,'vertical', 'LineStyle','none','Color', 'k');
% 
% % make fit
% % NHI
% [b_nh, m_nh, sigm_nh, sigb_nh, chi2, q] = fitexy(pltr_mydgl_nh_goodfiles, pltr_thissig_goodfiles, hi4pi_err_fields, mytoterr_pos)
% % fitter_nh = linear_fit(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles);
% 
% % apply fit
% % fit_x_nh = linspace(min(pltr_mydgl_nh_goodfiles),max(pltr_mydgl_nh_goodfiles));
% fit_x_nh = linspace(0,max(pltr_mydgl_nh_goodfiles));
% fit_y_nh = (m_nh*fit_x_nh + b_nh);
% 
% % plot fit
% fit_nh = plot(fit_x_nh,fit_y_nh,'k');
% 
% ylabel('EBL + DGL (nW/m^2/sr)')
% xlabel('Neutral Hydrogen Column Density [cm^{-2}]')
% figcnt = figcnt + 1;

% Plot correlation between 100m emission and galactic b
% figure(figcnt); clf
%
% % Scaling factor for 100 micron emission with galactic latitude (free parameter A not included b/c part of b(lambda)
% dl = (1 - 1.1 .* (0.61).*sqrt(sin(abs(myb(isgood)).*pi./180)));
%
% % p1 = scatter(abs(myb(isgood)),pltr_my100m_planck(isgood),'b'); % b vs 100m
% % p1 = scatter(sqrt(sin(abs(myb(isgood)*pi/180))),pltr_my100m_planck(isgood),'b'); % sqrt(sin(b)) vs 100m
% % p1 = scatter(abs(myb(isgood)),pltr_my100m_planck(isgood).*sqrt(sin(abs(myb(isgood)*pi/180))),'b'); % b vs 100m*sqrt(sin(b))
% p1 = scatter(abs(myb(isgood)),pltr_my100m_planck(isgood).*dl,'b'); % b vs 100m*dl
%
% hold on
% % i1 = scatter(abs(myb(isgood)),pltr_my100m_iris(isgood),'g'); % b vs 100m
% % i1 = scatter(sqrt(sin(abs(myb(isgood)*pi/180))),pltr_my100m_iris(isgood),'g'); % sqrt(sin(b)) vs 100m
% % i1 = scatter(abs(myb(isgood)),pltr_my100m_iris(isgood).*sqrt(sin(abs(myb(isgood)*pi/180))),'g'); % b vs 100m*sqrt(sin(b))
% i1 = scatter(abs(myb(isgood)),pltr_my100m_iris(isgood).*dl,'g'); % b vs 100m*dl
%
% % i2 = scatter(abs(myb(isgood)),pltr_my100m_iris_sfd(isgood),'m'); % b vs 100m
% % i2 = scatter(sqrt(sin(abs(myb(isgood)*pi/180))),pltr_my100m_iris_sfd(isgood),'m'); % sqrt(sin(b)) vs 100m
% % i2 = scatter(abs(myb(isgood)),pltr_my100m_iris_sfd(isgood).*sqrt(sin(abs(myb(isgood)*pi/180))),'m'); % b vs 100m*sqrt(sin(b))
% i2 = scatter(abs(myb(isgood)),pltr_my100m_iris_sfd(isgood).*dl,'m'); % b vs 100m*dl
%
% % idealRelation = sqrt(sin(abs(myb(isgood)*pi/180)));
% % idealRelation = (max(pltr_my100m_planck(isgood))-min(pltr_my100m_planck(isgood)))/(max(idealRelation)-min(idealRelation))*idealRelation + (max(pltr_my100m_planck(isgood)) - (max(pltr_my100m_planck(isgood))-min(pltr_my100m_planck(isgood)))/(max(idealRelation)-min(idealRelation))*max(idealRelation));
% % pT = plot(abs(myb(isgood)),idealRelation,'o');
% % i2 = scatter(abs(myb(isgood)),pltr_my100m_iris_sfd(isgood).*idealRelation,'m');
%
% % make fit
% % NHI
% % [b_nh, m_nh, sigm_nh, sigb_nh, chi2, q] = fitexy(pltr_mydgl_nh_goodfiles, pltr_thissig_goodfiles, hi4pi_err_fields, pltr_sem)
% fitter_iris_sfd = linear_fit(abs(myb(isgood)),pltr_my100m_iris_sfd(isgood).*dl);
% fitter_iris = linear_fit(abs(myb(isgood)),pltr_my100m_iris(isgood).*dl);
% fitter_planck = linear_fit(abs(myb(isgood)),pltr_my100m_planck(isgood).*dl);
%
% % apply fit
% fit_x = linspace(min(abs(myb(isgood))),max(abs(myb(isgood))));
% fit_y_iris_sfd = (fitter_iris_sfd.m*fit_x + fitter_iris_sfd.b);
% fit_y_iris = (fitter_iris.m*fit_x + fitter_iris.b);
% fit_y_planck = (fitter_planck.m*fit_x + fitter_planck.b);
%
% % plot fit
% plot(fit_x,fit_y_iris_sfd,'m');
% plot(fit_x,fit_y_iris,'g');
% plot(fit_x,fit_y_planck,'b');
%
% % ylabel('100 \mum (MJy/sr)')
% % ylabel('100 \mum (MJy/sr) * sin(b)^{1/2}')
% ylabel('100 \mum (MJy/sr) * d(b)')
%
% xlabel('Galactic b (deg)')
% legend([p1,i1,i2],{'Planck','Iris','Iris + SFD'},'Location','northeast');
% figcnt = figcnt + 1;

% %lil_100m <-4 of them
% pltr_my100m_planck_goodfiles.*dl
% pltr_my100m_iris_goodfiles.*dl
% pltr_my100m_iris_sfd_goodfiles.*dl
% pltr_mydgl_nh_goodfiles.*dl
% 
% %sig_x <-4 of them
% pltr_my100merr_combo_planck_goodfiles.*dl
% iris_err_fields.*dl
% iris_err_fields.*dl
% hi4pi_err_fields
% 
% pltr_thissig_goodfiles %lil_optf
% mytoterr_pos %sig_y
% myext_goodfiles %ext


% fontsize = 12;
% fontname = 'Calibri';
% fontweight = 'normal';
% figure(figcnt); clf
% tl = tiledlayout(5,4,'TileSpacing','tight','padding','tight');
% set(gcf,'Color',[1 1 1]) %make background pure white (good for editing)
% for jj = 1:length(goodfiles)
%     nexttile
% %     subplot(4,5,jj)
%     p1 = pcolor(field_snaps_data(:,:,jj));
%     set(p1, 'EdgeColor', 'none');
%     set(gca,'Color',[0 0.4470 0.7410]) %make nans matlab-blue
%     colormap(gca,'gray')
%     
%     % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
%     a = colorbar;
%     a.Label.String = 'Intensity [DN]';
%     pbaspect([1 1 1]);
%     xlabel('LORRI X Pixels','FontName', 'Calibri');
%     ylabel('LORRI Y Pixels');
%     caxis([-10,10]);
%     %         title(sprintf('%s',datastruct.header.rawfile));
%     % grid minor;
%     % title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
%     title(sprintf('Field %d',unique_fields(jj)));
%     set(gca,'YDir','normal');
% %     gca.Title.FontSize = fontsize;
% %     gca.Title.FontName = fontname;
% %     gca.Title.FontWeight = fontweight;
% %     gca.FontSize = fontsize;
% %     gca.FontName = fontname;
% %     gca.FontWeight = fontweight;
% %     gca.XLabel.FontSize = fontsize;
% %     gca.XLabel.FontName = fontname;
% %     gca.XLabel.FontWeight = fontweight;
% %     gca.YLabel.FontSize = fontsize;
% %     gca.YLabel.FontName = fontname;
% %     gca.YLabel.FontWeight = fontweight;
% %     gca.Colorbar.FontSize = fontsize;
% %     gca.Colorbar.FontName = fontname;
% %     gca.Colorbar.FontWeight = fontweight;
% %     gca.Colorbar.Label.FontSize = fontsize;
% %     gca.Colorbar.Label.FontName = fontname;
% %     gca.Colorbar.Label.FontWeight = fontweight;
% end
% figcnt = figcnt + 1;
% 
% figure(figcnt); clf
% tl = tiledlayout(5,4,'TileSpacing','tight','padding','tight');
% set(gcf,'Color',[1 1 1]) %make background pure white (good for editing)
% for jj = 1:length(goodfiles)
%     nexttile
% %     subplot(4,5,jj)
%     p1 = pcolor(field_snaps_cal(:,:,jj));
%     set(p1, 'EdgeColor', 'none');
%     set(gca,'Color',[0 0.4470 0.7410]) %make nans matlab-blue
%     colormap(gca,'gray')
%     
%     % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
%     a = colorbar;
%     a.Label.String = 'Intensity [nW m^{-2} sr^{-1}]';
%     pbaspect([1 1 1]);
%     xlabel('LORRI X Pixels');
%     ylabel('LORRI Y Pixels');
%     caxis([-150,150]);
%     %         title(sprintf('%s',datastruct.header.rawfile));
%     % grid minor;
%     % title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
%     title(sprintf('Field %d',unique_fields(jj)));
%     set(gca,'YDir','normal');
% %     gca.Title.FontSize = fontsize;
% %     gca.Title.FontName = fontname;
% %     gca.Title.FontWeight = fontweight;
% %     gca.FontSize = fontsize;
% %     gca.FontName = fontname;
% %     gca.FontWeight = fontweight;
% %     gca.XLabel.FontSize = fontsize;
% %     gca.XLabel.FontName = fontname;
% %     gca.XLabel.FontWeight = fontweight;
% %     gca.YLabel.FontSize = fontsize;
% %     gca.YLabel.FontName = fontname;
% %     gca.YLabel.FontWeight = fontweight;
% %     gca.Colorbar.FontSize = fontsize;
% %     gca.Colorbar.FontName = fontname;
% %     gca.Colorbar.FontWeight = fontweight;
% %     gca.Colorbar.Label.FontSize = fontsize;
% %     gca.Colorbar.Label.FontName = fontname;
% %     gca.Colorbar.Label.FontWeight = fontweight;
% end
% figcnt = figcnt + 1;

%---stacked bar plot or error constituents---
figure(figcnt); clf
h = bar([pltr_unc,myscattering_toterr_goodfiles,mypsfwing_psferr_goodfiles,...
    mytrierr_goodfiles,myghostdifferrpos_goodfiles,abs(mygalerr_goodfiles),abs(mymagerr_goodfiles),],'stacked');
h_barWidth = h.BarWidth;
h_barHeight = sum([pltr_unc,myscattering_toterr_goodfiles,mypsfwing_psferr_goodfiles,...
    mytrierr_goodfiles,myghostdifferrpos_goodfiles,abs(mygalerr_goodfiles),abs(mymagerr_goodfiles),],2);
for i=1:length(unique_fields)
    text(i, h_barHeight(i), num2str(i),'HorizontalAlignment','center','VerticalAlignment','bottom');
end 
ylabel('Cumulative Error [nW m^{-2} sr^{-1}]')
xlabel('Field Number')
legend('Statistical','Scattering','PSF Wing','TRILEGAL','Diffuse Ghosts','Masking Galaxies','Masking Stars',...
    'Location','bestoutside','Orientation','horizontal')
figcnt = figcnt + 1;

% figure(figcnt); clf
% h = bar([abs(mymagerr_goodfiles),abs(mygalerr_goodfiles),myghostdifferrpos_goodfiles,mytrierr_goodfiles,...
%     mypsfwing_psferr_goodfiles,myscattering_toterr_goodfiles,myunc],'stacked');
% h_barWidth = h.BarWidth;
% h_barHeight = sum([myunc,myscattering_toterr_goodfiles,mypsfwing_psferr_goodfiles,...
%     mytrierr_goodfiles,myghostdifferrpos_goodfiles,abs(mygalerr_goodfiles),abs(mymagerr_goodfiles),],2);
% for i=1:length(unique_fields)
%     text(i, h_barHeight(i), num2str(i),'HorizontalAlignment','center','VerticalAlignment','bottom');
% end 
% ylabel('Cumulative Error [nW m{^-2} sr^{-1}]')
% xlabel('Field Number')
% legend('Statistical','Scattering','PSF Wing','TRILEGAL','Diffuse Ghosts','Masking Galaxies','Masking Stars',...
%     'Location','bestoutside','Orientation','horizontal')
% set(gca,'YScale','log')
% figcnt = figcnt + 1;

%=============Planck=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_planck_goodfiles.*dl;
sig_x = pltr_my100merr_combo_planck_goodfiles.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),1);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'Planck');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability
planck_optp = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_optp');
planck_residual = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_residual');
planck_opt = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_opt');
planck_residual_prext = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_residual_prext');
planck_x = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_x');
planck_sigx = sig_x;
planck_sigy = mytoterr_pos;

%=============IRIS=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_iris_goodfiles.*dl;
sig_x = iris_err_fields.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),1);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'IRIS');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability
iris_optp = h5read([pwd,'/fit_results_data_',strrep('IRIS','/',''),'.h5'],'/lil_optp');
iris_residual = h5read([pwd,'/fit_results_data_',strrep('IRIS','/',''),'.h5'],'/lil_residual');
iris_x = h5read([pwd,'/fit_results_data_',strrep('IRIS','/',''),'.h5'],'/lil_x');

%=============IRIS/SFD=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_iris_sfd_goodfiles.*dl;
sig_x = iris_err_fields.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),1);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'IRIS/SFD');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability
iris_sfd_optp = h5read([pwd,'/fit_results_data_',strrep('IRIS/SFD','/',''),'.h5'],'/lil_optp');
iris_sfd_residual = h5read([pwd,'/fit_results_data_',strrep('IRIS/SFD','/',''),'.h5'],'/lil_residual');
iris_sfd_x = h5read([pwd,'/fit_results_data_',strrep('IRIS/SFD','/',''),'.h5'],'/lil_x');

%=============NHI=============
hi4pi_err = ((2.3e18)/5)/sqrt((1.13*(16.2)^2)/((17.4^2))); % 5-sig rms sensitivity (2.3e18 cm^-2) converted to 1-sig and from per hi4pi beam to per lorri image
hi4pi_err_fields = ones(length(goodfiles),1)*hi4pi_err;
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_mydgl_nh_goodfiles.*dl;
sig_x = hi4pi_err_fields.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),1);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'NHI');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability
nhi_optp = h5read([pwd,'/fit_results_data_',strrep('NHI','/',''),'.h5'],'/lil_optp');
nhi_residual = h5read([pwd,'/fit_results_data_',strrep('NHI','/',''),'.h5'],'/lil_residual');
nhi_x = h5read([pwd,'/fit_results_data_',strrep('NHI','/',''),'.h5'],'/lil_x');

%---cleanup---
% delete 'fit_info.txt' 'fit_info_txt.txt'

%%%%% Redo fit without d(b) dependence to get blam %%%%%

%=============Planck=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_planck_goodfiles;
sig_x = pltr_my100merr_combo_planck_goodfiles;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),2);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'Planck');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability
planck_optp_nodb = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_optp');
planck_residual_nodb = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_residual');
planck_opt_nodb = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_opt');
planck_residual_prext_nodb = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_residual_prext');
planck_x_nodb = h5read([pwd,'/fit_results_data_',strrep('Planck','/',''),'.h5'],'/lil_x');
planck_sigx_nodb = sig_x;
planck_sigy_nodb = mytoterr_pos;


%=============IRIS=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_iris_goodfiles;
sig_x = iris_err_fields;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),2);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'IRIS');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability

%=============IRIS/SFD=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_iris_sfd_goodfiles;
sig_x = iris_err_fields;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),2);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'IRIS/SFD');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability

%=============NHI=============
hi4pi_err = ((2.3e18)/5)/sqrt((1.13*(16.2)^2)/((17.4^2))); % 5-sig rms sensitivity (2.3e18 cm^-2) converted to 1-sig and from per hi4pi beam to per lorri image
hi4pi_err_fields = ones(length(goodfiles),1)*hi4pi_err;
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_mydgl_nh_goodfiles;
sig_x = hi4pi_err_fields;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i),2);
end
fclose(fid);
%---make text title file---
fid = fopen('fit_info_txt.txt','w');
fprintf(fid,'NHI');
fclose(fid);
%---run pthon---
% pyrunfile("fit_results.py")
system(['python ',pwd,'/fit_results.py']);
disp(' '); %space for readability

%---cleanup---
delete 'fit_info.txt' 'fit_info_txt.txt'

% IPD Plots
myipd_goodfiles = [0.12611579;...
0.33819679;...
0.32918186;...
0.18524416;...
0.35001141;...
0.36698896;...
0.076865947;...
0.33709999;...
0.10806821;...
0.28672479;...
0.32283905;...
0.16193704;...
0.0950076;...
0.084372334;...
0.099879038;...
0.12701996;...
0.10457724;...
0.05585325;...
0.062172902];

ipd_poserr = [1.2611579;...
3.3819679;...
3.2918186;...
1.8524416;...
3.5001141;...
3.6698896;...
0.76865947;...
3.3709999;...
1.0806821;...
2.8672479;...
3.2283905;...
1.6193704;...
0.950076;...
0.84372334;...
0.99879038;...
1.2701996;...
1.0457724;...
0.5585325;...
0.62172902];

ipd_negerr = [0.012611579;...
0.033819679;...
0.032918186;...
0.018524416;...
0.035001141;...
0.036698896;...
0.007686595;...
0.033709999;...
0.010806821;...
0.028672479;...
0.032283905;...
0.016193704;...
0.00950076;...
0.008437233;...
0.009987904;...
0.012701996;...
0.010457724;...
0.005585325;...
0.00621729];

figure(figcnt);
clf
errorbar(myipd_goodfiles,planck_optp,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
ylabel('Correlative COB [nW m^{-2} sr^{-1}]')
xlabel('IPD [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;

figure(figcnt);
clf
errorbar(myipd_goodfiles,planck_residual,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
ylabel('COB Residual [nW m^{-2} sr^{-1}]')
xlabel('IPD [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;

fitter = linear_fit(myipd_goodfiles,planck_residual,myunc_planck);
fitter.m % Get slope
sqrt(fitter.variancem) % Get slope err
corrcoef(myipd_goodfiles,planck_residual) % Get corr coef

figure(figcnt);
clf
errorbar(mydist,myipd_goodfiles,ipd_negerr,ipd_poserr,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
% [dist_sorted, dist_sortedIndexes] = sort(mydist);
% plot(dist_sorted,myipd_goodfiles(dist_sortedIndexes),'Color',[0.8500, 0.3250, 0.0980]);
% hold on;
% % plot(dist_sorted,myipd_goodfiles(dist_sortedIndexes)+ipd_poserr(dist_sortedIndexes),'Color','b');
% % plot(dist_sorted,myipd_goodfiles(dist_sortedIndexes)-ipd_negerr(dist_sortedIndexes),'Color','g');
% pp = patch([dist_sorted; flip(dist_sorted)],...
%     [myipd_goodfiles(dist_sortedIndexes)-ipd_negerr(dist_sortedIndexes);...
%     flip(myipd_goodfiles(dist_sortedIndexes)+ipd_poserr(dist_sortedIndexes))], [0.8500, 0.3250, 0.0980], 'EdgeColor', 'none');
% alpha(pp,.5)
ylabel('IPD [nW m^{-2} sr^{-1}]')
xlabel('Solar Distance [AU]')
set(gca,'YScale','log')
figcnt = figcnt + 1;

% COB vs. distance from Sun
figure(figcnt);
clf
errorbar(mysun_goodfiles,planck_optp,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
ylabel('Correlative COB [nW m^{-2} sr^{-1}]')
xlabel('Solar Distance')
figcnt = figcnt + 1;

% Exclusion time vs COB plot
all_times = [0;25;50;75;100;125;150;175;200;225;250;275;300;325;350;375;400];
cob = [22.34;21.88;22.26;21.66;21.64;21.63;21.98;21.29;19.63;19.86;19.09;18.89;19.94;19.86;20.34;20.72;19.1425]; % mean of COB from fit
cob_err = [1.435;1.435;1.215;1.255;1.255;1.2525;1.23;1.305;1.97;1.955;2.105;2.0025;2.2325;2.2425;2.1875;2.1325;2.3025]; % This is mean of sigma COB from fit
all_fields = [26;25;23;21;20;19;19;18;13;12;10;10;8;8;8;8;6];

figure(figcnt);
clf
errorbar(all_times,cob,cob_err,'.','MarkerSize',30,'MarkerEdge',[0.4940 0.1840 0.5560],'LineStyle','none','LineWidth',2,'Color',[0.4940 0.1840 0.5560]); %old error bars were myunc, based on std not sem
xlim([-5,405]);
ylabel('Correlative COB [nW m^{-2} sr^{-1}]')
yyaxis right
plot(all_times,all_fields,'color','#696969','LineWidth',2);
ax = gca;
ax.YAxis(2).Color = 'k';
ylabel('Number of Fields')
xlabel('Time Excluded from Sequence Start [s]')
figcnt = figcnt + 1;

%----- Settings for COB residual vs gal lat plots -----
FLG_override_residual = true; %overrides residual with calculated one
FLG_define_slope = true; %true uses user-defined slope value, false uses calc'd one
planck_slope_desired = 5.25; %blam 5.5 makes slope 0
planck_slope_desired_prext = 5.25; %blam 5.5 makes slope 0 (for pre-extinction), 5.25 for no d(b) - round to 5 since 5 better than 5.5
FLG_guaranteeFlat = false; %[too shoddy to use atm]only wroks if FLG_define_slope is false, uses a loop to guarantee calc'd slope/intercept makes the data flat
% COB residual vs gal lat
figure(figcnt); clf
if( FLG_override_residual )
    if( FLG_define_slope )
        %dictate a slope
        planck_residual_new = planck_optp - ( planck_slope_desired*planck_x + 23.31 );
        %remove mean
        planck_residual_new_mean = mean(planck_residual_new); %get mean as var to see
        planck_residual_new_intercept = 23.31 + planck_residual_new_mean;
        planck_residual_new = planck_residual_new - planck_residual_new_mean; %equivalent to 23.31-planck_residual_mean
    else
        %estimate ideal slope and intercept
        planck_optp_slope = polyfit(planck_x, planck_optp, 1);
%         [intercept, slope, ~, ~, ~, ~] = fitexy(planck_x, planck_optp, planck_sigx, planck_sigy );
        planck_residual_new = planck_optp - ( planck_optp_slope(1)*planck_x + planck_optp_slope(2) );
        planck_residual_new_intercept = planck_optp_slope(2);
        if( FLG_guaranteeFlat )
            dir = -1;
            fitter = linear_fit(abs(mygal),planck_residual_new,myunc_planck);
            fitter_prev = abs(fitter.m);
            while( fitter_prev > 1E-6)
                planck_residual_new = planck_optp - ( planck_optp_slope(1)*planck_x + planck_optp_slope(2) );
                fitter = linear_fit(abs(mygal),planck_residual_new,myunc_planck);
                if( fitter.m > fitter_prev )
                    dir = -dir; %switch direction
                end
                planck_optp_slope(1) = planck_optp_slope(1) + 1E-6*dir; %instead of static 1E-6 using fitter.m as a scaling factor
                fitter_prev = abs(fitter.m); %new prev for next loop
            end
            planck_residual_new_mean = mean(planck_residual_new); %get mean as var to see
            planck_residual_new_intercept = planck_residual_new_intercept + planck_residual_new_mean;
        end
        planck_slope_desired = planck_optp_slope(1);
    end
    %plot
    disp(['Ideal slope: ',num2str(planck_slope_desired),' | Ideal intercept: ',num2str(planck_residual_new_intercept)])
    me = errorbar(abs(mygal),planck_residual_new,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_new,myunc_planck);
    fitter.m % Get slope
    corrcoef(abs(mygal),planck_residual_new) % Get corr coef
else
    me = errorbar(abs(mygal),planck_residual,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual,myunc_planck);
    fitter.m % Get slope
    corrcoef(abs(mygal),planck_residual) % Get corr coef
end
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('COB Residual [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;


% No db COB residual vs gal lat
figure(figcnt); clf %blam 6 makes slope 0
if( FLG_override_residual )
    if( FLG_define_slope )
        %dictate a slope
        planck_slope_desired_nodb = planck_slope_desired;
        planck_residual_new_nodb = planck_optp_nodb - ( planck_slope_desired_nodb*planck_x_nodb + 22.57 );
        %remove mean
        planck_residual_new_mean_nodb = mean(planck_residual_new_nodb); %get mean as var to see
        planck_residual_new_intercept_nodb = 22.57 + planck_residual_new_mean_nodb;
        planck_residual_new_nodb = planck_residual_new_nodb - planck_residual_new_mean_nodb; %equivalent to 22.57-planck_residual_mean
    else
        %estimate ideal slope and intercept
        planck_optp_slope_nodb = polyfit(planck_x_nodb, planck_optp_nodb, 1);
%         [intercept, slope, ~, ~, ~, ~] = fitexy(planck_x_nodb, planck_optp_nodb, planck_sigx_nodb, planck_sigy_nodb );
        planck_residual_new_nodb = planck_optp_nodb - ( planck_optp_slope_nodb(1)*planck_x_nodb + planck_optp_slope_nodb(2) );
        planck_residual_new_intercept_nodb = planck_optp_slope_nodb(2);
        if( FLG_guaranteeFlat )
            dir = -1;
            fitter = linear_fit(abs(mygal),planck_residual_new_nodb,myunc_planck);
            fitter_prev = fitter.m;
            while(fitter_prev > 1E-6)
                planck_residual_new_nodb = planck_optp_nodb - ( planck_optp_slope_nodb(1)*planck_x_nodb + planck_optp_slope_nodb(2) );
                fitter = linear_fit(abs(mygal),planck_residual_new_nodb,myunc_planck);
                if( fitter.m > fitter_prev )
                    dir = -dir; %switch direction
                end
                planck_optp_slope_nodb(1) = planck_optp_slope_nodb(1) + abs(fitter.m)*dir; %instead of static 1E-6 using fitter.m as a scaling factor
                fitter_prev = fitter.m; %new prev for next loop
            end
            planck_residual_new_mean_nodb = mean(planck_residual_new_nodb); %get mean as var to see
            planck_residual_new_intercept_nodb = planck_residual_new_intercept_nodb(2) + planck_residual_new_mean_nodb;
        end
        planck_slope_desired_nodb = planck_optp_slope_nodb(1);
    end
    %plot
    disp(['[nodb] Ideal slope: ',num2str(planck_slope_desired_nodb),' | Ideal intercept: ',num2str(planck_residual_new_intercept_nodb)])
    me = errorbar(abs(mygal),planck_residual_new_nodb,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_new_nodb,myunc_planck);
    fitter.m % Get slope
else
    me = errorbar(abs(mygal),planck_residual_nodb,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_nodb,myunc_planck);
    fitter.m % Get slope
end
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('COB Residual no d(b) [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;

% COB residual PRE-EXTINCTION vs gal lat
figure(figcnt); clf
if( FLG_override_residual )
    if( FLG_define_slope )
        %dictate a slope
        planck_residual_prext_new = planck_opt - ( planck_slope_desired_prext*planck_x + 22.67 );
        %remove mean
        planck_residual_prext_new_mean = mean(planck_residual_prext_new); %get mean as var to see
        planck_residual_prext_new_intercept = 22.67 + planck_residual_prext_new_mean;
        planck_residual_prext_new = planck_residual_prext_new - planck_residual_prext_new_mean; %equivalent to 22.67-planck_residual_mean
    else
        %estimate ideal slope and intercept
        planck_opt_slope = polyfit(planck_x, planck_opt, 1);
        planck_residual_prext_new = planck_opt - ( planck_opt_slope(1)*planck_x + planck_opt_slope(2) );
        planck_residual_prext_new_intercept = planck_opt_slope(2);
        if( FLG_guaranteeFlat )
            dir = -1;
            fitter = linear_fit(abs(mygal),planck_residual_prext_new,myunc_planck);
            fitter_prev = fitter.m;
            while(fitter_prev > 1E-6)
                planck_residual_prext_new = planck_opt - ( planck_opt_slope(1)*planck_x + planck_opt_slope(2) );
                fitter = linear_fit(abs(mygal),planck_residual_prext_new,myunc_planck);
                if( fitter.m > fitter_prev )
                    dir = -dir; %switch direction
                end
                planck_opt_slope(1) = planck_opt_slope(1) + abs(fitter.m)*dir; %instead of static 1E-6 using fitter.m as a scaling factor
                fitter_prev = fitter.m; %new prev for next loop
            end
            planck_residual_prext_new_mean = mean(planck_residual_prext_new); %get mean as var to see
            planck_residual_prext_new_intercept = planck_residual_prext_new_intercept(2) + planck_residual_prext_new_mean;
        end
        planck_slope_desired_prext = planck_opt_slope(1);
    end
    %plot
    disp(['[Pre-ext] Ideal slope: ',num2str(planck_slope_desired_prext),' | Ideal intercept: ',num2str(planck_residual_prext_new_intercept)])
    me = errorbar(abs(mygal),planck_residual_prext_new,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.4500, 0.3250, 0.2980],'LineStyle','none','Color',[0.4500, 0.3250, 0.2980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_prext_new,myunc_planck);
    fitter.m % Get slope
    corrcoef(abs(mygal),planck_residual_prext_new) % Get corr coef    
else
    me = errorbar(abs(mygal),planck_residual_prext,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.4500, 0.3250, 0.2980],'LineStyle','none','Color',[0.4500, 0.3250, 0.2980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_prext,myunc_planck);
    fitter.m % Get slope
    corrcoef(abs(mygal),planck_residual_prext) % Get corr coef
end
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('COB Residual Pre-Extinction [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;

% No db COB residual PRE-EXTINCTION vs gal lat
figure(figcnt); clf %blam 6 makes slope 0
% planck_residual_prext_nodb = planck_opt_nodb - ( planck_prext_slope_desired*planck_x_nodb + 22.15 );
if( FLG_override_residual )
    if( FLG_define_slope )
        %dictate a slope
        planck_slope_desired_prext_nodb = planck_slope_desired_prext;
        planck_residual_prext_new_nodb = planck_opt_nodb - ( planck_slope_desired_prext_nodb*planck_x_nodb + 22.15 );
        %remove mean
        planck_residual_prext_new_mean_nodb = mean(planck_residual_prext_new_nodb); %get mean as var to see
        planck_residual_prext_new_intercept_nodb = 22.15 + planck_residual_prext_new_mean_nodb;
        planck_residual_prext_new_nodb = planck_residual_prext_new_nodb - planck_residual_prext_new_mean_nodb; %equivalent to 22.15-planck_residual_mean
    else
        %estimate ideal slope and intercept
        planck_opt_slope_nodb = polyfit(planck_x_nodb, planck_opt_nodb, 1);
%         planck_opt_slope_nodb
%         [intercept, slope, ~, ~, ~, ~] = fitexy(planck_x_nodb, planck_opt_nodb, planck_sigx_nodb, planck_sigy_nodb );
%         [slope , intercept]
        planck_residual_prext_new_nodb = planck_opt_nodb - ( planck_opt_slope_nodb(1)*planck_x_nodb + planck_opt_slope_nodb(2) );
        planck_residual_prext_new_intercept_nodb = planck_opt_slope_nodb(2);
        if( FLG_guaranteeFlat )
            dir = -1;
            fitter = linear_fit(abs(mygal),planck_residual_prext_new_nodb,myunc_planck);
            fitter_prev = fitter.m;
            while(fitter_prev > 1E-6)
                planck_residual_prext_new_nodb = planck_opt_nodb - ( planck_opt_slope_nodb(1)*planck_x_nodb + planck_opt_slope_nodb(2) );
                fitter = linear_fit(abs(mygal),planck_residual_prext_new_nodb,myunc_planck);
                if( fitter.m > fitter_prev )
                    dir = -dir; %switch direction
                end
                planck_opt_slope_nodb(1) = planck_opt_slope_nodb(1) + abs(fitter.m)*dir; %instead of static 1E-6 using fitter.m as a scaling factor
                fitter_prev = fitter.m; %new prev for next loop
            end
            planck_residual_prext_new_mean_nodb = mean(planck_residual_prext_new_nodb); %get mean as var to see
            planck_residual_prext_new_intercept_nodb = planck_residual_prext_new_intercept_nodb(2) + planck_residual_prext_new_mean_nodb;
        end
        planck_slope_desired_prext_nodb = planck_opt_slope_nodb(1);
    end
    %plot
    disp(['[nodb pre-ext] Ideal slope: ',num2str(planck_slope_desired_prext_nodb),' | Ideal intercept: ',num2str(planck_residual_prext_new_intercept_nodb)])
    me = errorbar(abs(mygal),planck_residual_prext_new_nodb,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.4500, 0.3250, 0.2980],'LineStyle','none','Color',[0.4500, 0.3250, 0.2980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_prext_new_nodb,myunc_planck);
    fitter.m % Get slope
else
    me = errorbar(abs(mygal),planck_residual_prext_nodb,myunc_planck,'.','MarkerSize',20,'MarkerEdge',[0.4500, 0.3250, 0.2980],'LineStyle','none','Color',[0.4500, 0.3250, 0.2980]); %old error bars were myunc, based on std not sem
    fitter = linear_fit(abs(mygal),planck_residual_prext_nodb,myunc_planck);
    fitter.m % Get slope
end
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('COB Residual Pre-Extinction no d(b) [nW m^{-2} sr^{-1}]')
figcnt = figcnt + 1;



% Planck 100m vs gal lat
figure(figcnt); clf 
me = errorbar(abs(mygal),pltr_my100m_planck_goodfiles,pltr_my100merr_combo_planck_goodfiles,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('Planck 100m Intensity [MJy sr^{-1}]')
figcnt = figcnt + 1;

% Planck 100m * db vs gal lat
figure(figcnt); clf 
me = errorbar(abs(mygal),pltr_my100m_planck_goodfiles.*dl,pltr_my100merr_combo_planck_goodfiles.*dl,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
xlabel(['Galactic Latitude [' char(176) ']'])
ylabel('Planck 100m Intensity [MJy sr^{-1}]')
figcnt = figcnt + 1;

distance = mydist;
rawmean = mysubmen;
rawerr = mysuberr;
cobmean = mymean;
coberr = myunc;
supermean = supermean;
supererr = superunc;
ohm = myohmp;

save('../scratch/nh_make_results.mat','distance','rawmean','rawmean',...
    'rawerr','cobmean','coberr','supermean','supererr','ohm')

% end

function parsave_image(fname, image) %this lets parfor save
  save(fname, 'image')
end
function parsave_nanimage(fname, nanimage) %this lets parfor save
  save(fname, 'nanimage')
end
