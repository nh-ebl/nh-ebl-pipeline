function nh_make_resultsCombo()

dglparams = nh_get_dgl_params();

want_errmags = 0;
want_errpsf = 0;

%Import paths for data location
paths = get_paths_old();
npaths = get_paths_new();
lpaths = get_paths_lauer();
wpaths = get_paths_newest();

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
oldgoodfields = [3,5,6,7];
% oldgoodfields = []; % Set good fields to none to skip data set
newgoodfields = [5,6,7,8];
newgood_exlude_enable = false; %enables skipping of sequences
% newgoodfields = [];
lauergoodfields = [1,2,3,4,5,6,7];
lauer_exlude_enable = false; %enables skipping of new sequences
% lauergoodfields = [];
newestgoodfields = [2,4,5,6,7,12,15,16,17,19,20,23];
newest_exlude_enable = false; %enables skipping of new sequences
% newestgoodfields = [];


% Calculate total fields
goodfiles = [oldgoodfields,newgoodfields,lauergoodfields,newestgoodfields];
zones = cumsum([1,length(oldgoodfields),length(newgoodfields),length(lauergoodfields),length(newestgoodfields)]); %Record the seperate zones and their orders
unique_fields = [1:length(goodfiles)];

%Check for old light files
parpoolobj = gcp('nocreate'); % check for thread pool, which can't use the load call
if isa(parpoolobj,'parallel.ThreadPool')
    delete(parpoolobj); %threads can't use load and will error out
end
parfor ifile=1:numel(lightfiles)
    datatemp = load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == oldgoodfields)
        isgoodold(ifile) = 1;
    end
end

%Check for new light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = 150; %s, time to skip at start of sequence
newgood_exclude = zeros(numel(llightfiles),1); %will fill up with new sequences to ignore
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
end
newgood_exclude = find(newgood_exclude); %get it into indexes

%Check for lauer light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = 150; %s, time to skip at start of sequence
lauer_exclude = zeros(numel(llightfiles),1); %will fill up with new sequences to ignore
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
end
lauer_exclude = find(lauer_exclude); %get it into indexes

%Check for newest light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = 150; %s, time to skip at start of sequence
newest_exclude = zeros(numel(wlightfiles),1); %will fill up with new sequences to ignore
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
end
% disp(['Total excluded: ',num2str(sum(newest_exclude))])
newest_goodfilesNum = sum(isgoodnewest&newest_exclude);
disp(['Newest|Good over files after isgoodnewest (good field check) and newest_exclude (time skip): ',num2str(newest_goodfilesNum)])
newest_exclude = find(newest_exclude); %get it into indexes

load('run_params.mat','params')

% If method file exists, read saved text file for data_type and see which method last used
flag_method = 'new';

%Number of old and new light files corresponding to good fields
numoldlightfiles = sum(isgoodold);
numnewlightfiles = sum(isgoodnew);
numlauerlightfiles = sum(isgoodlauer);
numnewestlightfiles = sum(isgoodnewest);

%Preallocate
isgood = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
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

%For old data files
fprintf('Loading old data \n');
jfile = 1;
% k = 1;
zones_indexer = 1; %Keeps count of the current zone
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
for ifile=1:numel(lightfiles)
    %If file is for a good field, load and save values
    if isgoodold(ifile) == 1
        %Load data files
        load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(lightfiles,1)));

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

        myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);

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
        pltr_my100merr_iris(jfile) = data.dgl.ohmstd_iris;
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
                myscattering_toterr(jfile) = sqrt(myscattering_gaiaerr(jfile)^2 + myscattering_masanaerr(jfile)^2);
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
                save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
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
            save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
        end
        jfile = jfile + 1;
    end
end

%For new data files
fprintf('Loading new data \n');
jfile = 1+numoldlightfiles; %redundant but kept to deal with exclusions
% k = k + 1;
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
for ifile=1:numel(nlightfiles)
    %If file is for a good field, load and save values
    if (isgoodnew(ifile) == 1) && ~any(newgood_exclude == ifile)

        %Load data files
        load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(nlightfiles,1)));

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

        myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);

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
        pltr_my100merr_iris(jfile) = data.dgl.ohmstd_iris;
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
                myscattering_toterr(jfile) = sqrt(myscattering_gaiaerr(jfile)^2 + myscattering_masanaerr(jfile)^2);
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
                save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
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
            save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
        end
        jfile = jfile + 1;
    end
end

%For lauer data files
fprintf('Loading Lauer data \n');
jfile = 1+numoldlightfiles+numnewlightfiles; %redundant but kept in case exclude is on
% k = k + 1;
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
for ifile=1:numel(llightfiles)
    %If file is for a good field and not excluded, load and save values
    if (isgoodlauer(ifile) == 1) && ~any(lauer_exclude == ifile)

        %Load data files
        load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(llightfiles,1)));

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

        myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);

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
        pltr_my100merr_iris(jfile) = data.dgl.ohmstd_iris;
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
                myscattering_toterr(jfile) = sqrt(myscattering_gaiaerr(jfile)^2 + myscattering_masanaerr(jfile)^2);
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
                save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
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
            save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
        end
        jfile = jfile + 1;
    end
end

%For newest data files
fprintf('Loading newest data \n');
jfile = 1+numoldlightfiles+numnewlightfiles+numlauerlightfiles; %redundant but kept in case exclude is on
% k = k + 1;
zones_indexer = zones_indexer + 1; %Increment the zone to the next zone of data
zones_goodfiles = goodfiles(zones(zones_indexer):zones(zones_indexer+1)-1); %Current zone's good fields
for ifile=1:numel(wlightfiles)
    %If file is for a good field and not excluded, load and save values
    if (isgoodnewest(ifile) == 1) && ~any(newest_exclude == ifile)

        %Load data files
        load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
        disp(sprintf('On file %d of %d.',ifile,size(wlightfiles,1)));

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

        myisl(jfile) = data.isl.trimean  + mypsfwing(jfile);

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
        pltr_my100merr_iris(jfile) = data.dgl.ohmstd_iris;
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
                myscattering_toterr(jfile) = sqrt(myscattering_gaiaerr(jfile)^2 + myscattering_masanaerr(jfile)^2);
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
                save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
            else
                newest_goodfilesNum = newest_goodfilesNum - 1; %reduce by 1 if there was a bad file
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
            save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(jfile),jfile),'image');

            nanimage = image;
            nanimage(data.mask.onemask) = NaN;
            save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(jfile),jfile),'nanimage');
        end
        jfile = jfile + 1;
    end
end
disp(['Newest|Total good files to be used (including good field (isgoodnewest), time limit (newest_exclude), isgood (header.bad check)): ',num2str(newest_goodfilesNum)])
% Select only data marked as good
isgood = logical(isgood);

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
    if any(k_field)
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

figure(1); clf
plot(mysun(isgood),mysig(isgood)-myisl(isgood)-pltr_mydgl_planck(isgood),'o')
xlabel('Solar Distance')
ylabel('EBL')

if strcmp(flag_method,'new') == 1
    figure(2); clf
    plot(mysun(isgood),myghostdiff(isgood),'o')
%     hold on;
%     errorbar(mysun(isgood),myghostdiff(isgood), myghostdifferrneg(isgood), myghostdifferrpos(isgood), 'LineStyle','none','Color','k');
    xlabel('Solar Distance')
    ylabel('Summed Ghost Intensity [nW m^{-2} sr^{-1}]')

    %     figure(3); clf
    %     plot(mysun(isgood),myghostdiffcomp(isgood),'o')
    %     yline(mean(myghostdiffcomp(isgood)))
    %     xlabel('Solar Distance')
    %     ylabel('(Predicted Summed Ghost Intensity - ISL Estimation) [nW m^{-2} sr^{-1}]')
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

% Calculate per-field values
for ifield=1:numel(goodfiles)

    whpl = unique_field_id == unique_fields(ifield);
    mydist(ifield) = mean(mysun(whpl & isgood));
    mygal(ifield) = mean(myb(whpl & isgood));
    if strcmp(flag_method,'new') == 1
        thissig = mysig(whpl & isgood)-myisl(whpl & isgood)-pltr_mydgl_planck(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        pltr_totsig_goodfiles = mysig(whpl & isgood)-myisl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        pltr_thissig_goodfiles(ifield) = mean(mysig(whpl & isgood)-myisl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood));
    else
        thissig = mysig(whpl & isgood)-myisl(whpl & isgood)-pltr_mydgl_planck(whpl & isgood);
    end
    mysig_goodfiles(ifield) = mean(mysig(whpl & isgood));
    mysig_magerr_goodfiles(ifield) = mean(mysig_magerr(whpl & isgood));
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
    mymean(ifield) = sum(thissig./thiserr.^2)./sum(1./thiserr.^2);
    myunc(ifield) = std(thissig);
    mysem(ifield) = std(thissig)/sqrt(length(thissig));
    pltr_sem(ifield) = std(pltr_totsig_goodfiles)/sqrt(length(pltr_totsig_goodfiles));
    %     myavp(ifield) = sum(myav(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    %     myohmp(ifield) = sum(myohm(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    pltr_mydataset(ifield) = median(mydataset(whpl)); %for plotting per dataset color
%     mytoterr_pos(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
%         mytrierr_goodfiles(ifield)^2 + myghostdifferrpos_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + mysem(ifield)^2);
%     mytoterr_neg(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
%         mytrierr_goodfiles(ifield)^2 + myghostdifferrneg_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + mysem(ifield)^2);
    
    % Use std (myunc) instead of sem (mysem) for individual statistical
    % errors
    mytoterr_pos(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
        mytrierr_goodfiles(ifield)^2 + myghostdifferrpos_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + myunc(ifield)^2);
    mytoterr_neg(ifield) = sqrt(mygalerr_goodfiles(ifield)^2 + mymagerr_goodfiles(ifield)^2 + mypsferr_goodfiles(ifield)^2 + ...
        mytrierr_goodfiles(ifield)^2 + myghostdifferrneg_goodfiles(ifield)^2 + myscattering_toterr_goodfiles(ifield)^2 + myunc(ifield)^2);
end
pltr_mydatasetTriplets(pltr_mydataset == 1,:) = repmat([0, 0.4470, 0.7410],sum(pltr_mydataset == 1),1); %set consistent color triplets
pltr_mydatasetTriplets(pltr_mydataset == 2,:) = repmat([0.8500, 0.3250, 0.0980],sum(pltr_mydataset == 2),1);
pltr_mydatasetTriplets(pltr_mydataset == 3,:) = repmat([0.9290, 0.6940, 0.1250],sum(pltr_mydataset == 3),1);

figure(4); clf
plot(mydist,mysubmen,'o')
hold on
errorbar(mydist,mysubmen,mysuberr,'o');
xlabel('Solar Distance')
ylabel('Error-Weighted Image Mean')

mysubmen

figure(5); clf
FLG_perDatasetColoring = true; %false plots default (same color), true plots colors per dataset (old/new/lar)
%plot(mydist,mymean,'.','MarkerSize',0);
%hold on
if ~FLG_perDatasetColoring
    errorbar(mydist,mymean,mysem,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
else
    hold on;
    for j = 1:length(mydist)
        errorbar(mydist(j),mymean(j),mysem(j),'.','MarkerSize',20,'MarkerEdge',pltr_mydatasetTriplets(j,:),'LineStyle','none','Color',pltr_mydatasetTriplets(j,:)); %old error bars were myunc, based on std not sem
    end
end
hold on;
xlabel('Solar Distance')
ylabel('Error-Weighted EBL')

supermean = sum(mymean./mysem.^2)./sum(1./mysem.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./mysem.^2)) %old way was myunc instead of mysem

% superav = sum(myavp./myunc.^2)./sum(1./myunc.^2)

supersub = sum(mysubmen./myunc.^2)./sum(1./myunc.^2)

plot(xlim,[supermean,supermean],'Color',[0, 0.4470, 0.7410])
plot(xlim,[supermean+superunc,supermean+superunc],'Color',[0, 0.4470, 0.7410],'LineStyle',':')
plot(xlim,[supermean-superunc,supermean-superunc],'Color',[0, 0.4470, 0.7410],'LineStyle',':')
plot(xlim,[0,0],'k:')
xlabel('Heliocentric Distance (AU)')
ylabel('Optical EBL (nW/m^2/sr)')
%   title('Zemcov et al, Preliminary w/ Stat. Errors')

figure(6); clf
me = errorbar(mygal,mymean,mysem,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
hold on;

z_pts = [15.4,18.1,-6.3,29.2];
z_err = [14.4,26.2,9.1,20.5];
z_err_sem = [14.4/sqrt(10),26.2/sqrt(10),9.1/sqrt(3),20.5/sqrt(3)];
z_b = [85.74,28.41,57.69,62.03];

zem = errorbar(z_b,z_pts,z_err_sem,'.','MarkerSize',20,'MarkerEdge',[0, 0.4470, 0.7410],'LineStyle','none','Color',[0, 0.4470, 0.7410]); %old error bars were myunc, based on std not sem
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
xlabel(['Galactic Latitude (' char(176) ')'])
ylabel('Optical EBL (nW/m^2/sr)')

% Plot DGL vs. EBL + DGL for all images
figure(7); clf
p1 = scatter(pltr_mydgl_planck(isgood),pltr_thissig,'o');
hold on
% p1err = errorbar(pltr_thissig,pltr_mydgl_planck,pltr_mydglerr_planck,'o');
i1 = scatter(pltr_mydgl_iris(isgood),pltr_thissig,'g');
% i1err = errorbar(pltr_thissig,pltr_mydgl_iris,pltr_mydglerr_iris,'m');
ylabel('EBL + DGL (nW/m^2/sr)')
xlabel('DGL (nW/m^2/sr)')
legend([p1,i1],{'Planck','Iris'},'Location','northwest');

% Plot DGL vs. EBL + DGL per field
figure(8); clf
p1 = scatter(pltr_mydgl_planck_goodfiles,pltr_thissig_goodfiles,'o');
hold on
% p1err = errorbar(pltr_thissig,pltr_mydgl_planck,pltr_mydglerr_planck,'o');
i1 = scatter(pltr_mydgl_iris_goodfiles,pltr_thissig_goodfiles,'g');
% i1err = errorbar(pltr_thissig,pltr_mydgl_iris,pltr_mydglerr_iris,'m');
ylabel('EBL + DGL (nW/m^2/sr)')
xlabel('DGL (nW/m^2/sr)')
legend([p1,i1],{'Planck','Iris'},'Location','northwest');

% Plot 100 micron emission vs. EBL + DGL with fit that includes x and y errors
figure(9); clf
FLG_perDatasetColoring = false; %false plots colors per Planck/Iris/Iris+SDF, true plots colors per dataset (old/new/lar)

% IRIS error
iris_err = (0.06)/sqrt((1.13*(4.3)^2)/((17.4^2))); % rms noise (0.06 MJy/sr) converted from per iris beam to per lorri image
iris_err_fields = ones(length(goodfiles),1)*iris_err;
% Scaling factor for 100 micron emission with galactic latitude (free parameter A not included b/c part of b(lambda)
dl = (1 - 1.1 .* (dglparams.g(1)).*sqrt(sin(abs(myb_goodfiles).*pi./180)));
% Calculate propagated error on d(b) due to error on g
dberr = (dglparams.g(2).^2 .* (1.1 .*sqrt(sin(abs(myb_goodfiles).*pi./180))).^2);
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

% Plot NHI vs. EBL + DGL
% figure(10); clf
%
% % HI4PI error
hi4pi_err = ((2.3e18)/5)/sqrt((1.13*(16.2)^2)/((17.4^2))); % 5-sig rms sensitivity (2.3e18 cm^-2) converted to 1-sig and from per hi4pi beam to per lorri image
hi4pi_err_fields = ones(length(goodfiles),1)*hi4pi_err;
% n1 = scatter(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,'k');
% hold on
% n1errx = errorbar(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,hi4pi_err_fields,'horizontal', 'LineStyle','none','Color', 'k');
% n1erry = errorbar(pltr_mydgl_nh_goodfiles,pltr_thisiris_err_fields.*dlsig_goodfiles,pltr_sem,'vertical', 'LineStyle','none','Color', 'k');
%
% % make fit
% % NHI
% [b_nh, m_nh, sigm_nh, sigb_nh, chi2, q] = fitexy(pltr_mydgl_nh_goodfiles, pltr_thissig_goodfiles, hi4pi_err_fields, pltr_sem)
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
% xlabel('Neutral Hydrogen Column Density (cm^{-2})')

% Plot correlation between 100m emission and galactic b
% figure(11); clf
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

%---stacked bar plot or error constituents---
figure(12); clf
bar([myunc,myscattering_toterr_goodfiles,mypsfwing_psferr_goodfiles,...
    mytrierr_goodfiles,myghostdifferrpos_goodfiles,mygalerr_goodfiles,abs(mymagerr_goodfiles),],'stacked')
ylabel('Cumulative Error [nW m{^-2} sr^{-1}]')
xlabel('Field Number')
legend('Statistical','Scattering','PSF Wing','TRILEGAL','Diffuse Ghosts','Masking Galaxies','Masking Stars',...
    'Location','bestoutside','Orientation','horizontal')


figure(13); clf
bar([abs(mymagerr_goodfiles),mygalerr_goodfiles,myghostdifferrpos_goodfiles,mytrierr_goodfiles,...
    mypsfwing_psferr_goodfiles,myscattering_toterr_goodfiles,myunc],'stacked')
ylabel('Cumulative Error [nW m{^-2} sr^{-1}]')
xlabel('Field Number')
legend('Statistical','Scattering','PSF Wing','TRILEGAL','Diffuse Ghosts','Masking Galaxies','Masking Stars',...
    'Location','bestoutside','Orientation','horizontal')
set(gca,'YScale','log')

%=============Planck=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_planck_goodfiles.*dl;
sig_x = pltr_my100merr_combo_planck_goodfiles.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i));
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

%=============IRIS=============
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_my100m_iris_goodfiles.*dl;
sig_x = iris_err_fields.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i));
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
lil_100m = pltr_my100m_iris_sfd_goodfiles.*dl;
sig_x = iris_err_fields.*dl;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i));
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
%---make data file---
fid = fopen('fit_info.txt','w');
lil_100m = pltr_mydgl_nh_goodfiles.*dl;
sig_x = hi4pi_err_fields;
for i=1:length(lil_100m)
    %Order is lil_opt, lil_100m, sig_y, sig_x, ext [5 total]
   fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',pltr_thissig_goodfiles(i),lil_100m(i),mytoterr_pos(i),sig_x(i),myext_goodfiles(i));
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
