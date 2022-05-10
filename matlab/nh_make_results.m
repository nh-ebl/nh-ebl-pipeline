function nh_make_results()
dglparams = nh_get_dgl_params();

want_errmags = 0;
want_errpsf = 0;
want_errtri = 0;
want_errdgl = 0;

% Set desired source of DGL calculation
% dgl_type = 'planck';
dgl_type = 'iris';
% dgl_type = 'iris_sfd';

% paths = get_paths_old();
% paths = get_paths_new();
paths = get_paths_lauer();

load('run_params.mat','params')

% Save which data we're looking at
if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/good/') == 1
    data_type = 'old';
    old = 1;
    new = 0;
elseif strcmp(paths.datadir,'/data/symons/nh_data/mat/') == 1
    data_type = 'new';
    old = 0;
    new = 1;
elseif strcmp(paths.datadir,'/data/symons/nh_data_lauer/mat/') == 1
    data_type = 'new';
    old = 0;
    new = 0;
    lauer = 1;
end

% If method file exists, read saved text file for data_type and see which method last used
fileID = fopen([data_type,'method.txt'],'r');
flag_method = fscanf(fileID,'%s');
fclose(fileID);

datafiles = dir(sprintf('%s*.mat',paths.datadir));

nfiles = size(datafiles,1);

isgood = zeros(nfiles,1);
mydate = zeros(nfiles,1);
myfieldnum = zeros(nfiles,1);
mytarget = cellstr('');
mysun = zeros(nfiles,1);
mydt = zeros(nfiles,1);
myelong = zeros(nfiles,1);
myexp = zeros(nfiles,1);
myref = zeros(nfiles,1);
myeng = zeros(nfiles,1);
myohm = zeros(nfiles,1);
myav = zeros(nfiles,1);
mysig = zeros(nfiles,1);
mymasked = zeros(nfiles,1);
biaslevl = zeros(nfiles,1);
biasmthd = zeros(nfiles,1);
biasdiff = zeros(nfiles,1);
rsolar = zeros(nfiles,1);
myb = zeros(nfiles,1);
mytri = zeros(nfiles,1);
mytrierr = zeros(nfiles,1);
mypsfwing = zeros(nfiles,1);
myghost = zeros(nfiles,1);
myscattering_tot = zeros(nfiles,1);
myisl = zeros(nfiles,1);
mydgl = zeros(nfiles,1);
mydglerr = zeros(nfiles,1);
myerr = zeros(nfiles,1);
mycrr = zeros(nfiles,1);
mymaskmean = zeros(nfiles,1);
myext = zeros(nfiles,1);
myghostdiff = zeros(nfiles,1);
myghostdifferrpos = zeros(nfiles,1);
myghostdifferrneg = zeros(nfiles,1);
myghostdiffcomp = zeros(nfiles,1);
pltr_mydgl_planck = zeros(nfiles,1);
pltr_mydglerr_planck = zeros(nfiles,1);
% pltr_myohmim_planck = zeros(nfiles,1);
pltr_mydgl_iris = zeros(nfiles,1);
pltr_mydglerr_iris = zeros(nfiles,1);
% pltr_myohmim_iris = zeros(nfiles,1);
pltr_mydgl_iris_sfd = zeros(nfiles,1);
pltr_mydglerr_iris_sfd = zeros(nfiles,1);
pltr_my100m_planck = zeros(nfiles,1);
pltr_my100merr_planck = zeros(nfiles,1);
pltr_my100merr_beta_planck = zeros(nfiles,1);
pltr_my100merr_temp_planck = zeros(nfiles,1);
pltr_my100merr_tau_planck = zeros(nfiles,1);
pltr_my100merr_combo_planck = zeros(nfiles,1);
pltr_my100m_iris = zeros(nfiles,1);
pltr_my100merr_iris = zeros(nfiles,1);
pltr_my100m_iris_sfd = zeros(nfiles,1);
pltr_my100merr_iris_sfd = zeros(nfiles,1);
pltr_mydgl_nh = zeros(nfiles,1);
pltr_mydglerr_nh = zeros(nfiles,1);
% set good files based on data set
if new == 1
    goodfiles = [5,6,7,8];
elseif old == 1
    goodfiles = [3,5,6,7];
elseif lauer == 1
    goodfiles = [1,2,3,4,5,6,7];
end

numfields = length(goodfiles);

thissum = 0;

for ifile=1:nfiles
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    mydate(ifile) = data.header.date_jd-data.header.launch_jd;
    
    myfieldnum(ifile) = data.header.fieldnum;
    
    mytarget{ifile} = data.header.date_cal;
    
    mysun(ifile) = data.distance.d_sun_NH;
    
    mydt(ifile) = data.distance.d_NH_target;%./1.496e8;;
    
    myelong(ifile) = data.coords.sol_elon;
    
    myb(ifile) = data.coords.galactic(2);
    
    myexp(ifile) = round(data.header.exptime);
    
    myref(ifile) = data.ref.mean;
    
    myeng(ifile) = data.ref.engmean;
    
%     myohm(ifile) = data.dgl.ohmmean;
    
    myav(ifile) = 0.4.*data.header.A_V;
    
    myext(ifile) = data.ext.mean_flux_rat;
    
    % Reference-corrected mean of masked image
    if want_errmags == 1
        mysig(ifile) = data.err_mags.stats.corrmean;
    else
        mysig(ifile) = data.stats.corrmean;
    end
    
    mymaskmean(ifile) = data.stats.maskmean;
    
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
    rsolar(ifile) = data.astrom.rsolar;

    % Reference-corrected most probable value of masked image
    %     mysig(ifile) = data.stats.corrmostprob;
    
    % Sum of image brightness of *masked* stars - how bright is the
    % masked-out portion of each image per pixel
    if want_errmags == 1
        starmaskfrac = sum(data.err_mags.mask.starmask(:))./256.^2;
        mymasked(ifile) = (sum(sum(data.err_mags.image.calimage.*data.err_mags.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
    else
        starmaskfrac = sum(data.mask.starmask(:))./256.^2;
        mymasked(ifile) = (sum(sum(data.image.calimage.*data.mask.starmask)))/(data.astrom.imagew.*data.astrom.imageh.*starmaskfrac);
    end
    
    mytri(ifile) = data.isl.trimean;
    mytrierr(ifile) = data.isl.trierr;
    
    if want_errpsf == 1
        mypsfwing(ifile) = data.err_psf.isl.usnowing;
    else
        mypsfwing(ifile) = data.isl.usnowing;
    end
    
    myisl(ifile) = data.isl.trimean  + mypsfwing(ifile);
    
    if strcmp(dgl_type,'planck') == 1
        mydgl(ifile) = data.dgl.dglmean_planck;
        mydglerr(ifile) = data.dgl.dglerr_planck;
    elseif strcmp(dgl_type, 'iris') == 1
        mydgl(ifile) = data.dgl.dglmean_iris;
        mydglerr(ifile) = data.dgl.dglerr_iris;
    elseif strcmp(dgl_type, 'iris_sfd') == 1
        mydgl(ifile) = data.dgl.dglmean_iris_sfd;
        mydglerr(ifile) = data.dgl.dglerr_iris_sfd;
    end
    pltr_mydgl_planck(ifile) = data.dgl.dglmean_planck;
    pltr_mydglerr_planck(ifile) = data.dgl.dglerr_planck;
%     pltr_myohmim_planck(ifile) = data.dgl.ohmim_planck;
    pltr_mydgl_iris(ifile) = data.dgl.dglmean_iris;
    pltr_mydglerr_iris(ifile) = data.dgl.dglerr_iris;
%     pltr_myohmim_iris(ifile) = data.dgl.ohmim_iris;
    pltr_mydgl_iris_sfd(ifile) = data.dgl.dglmean_iris_sfd;
    pltr_mydglerr_iris_sfd(ifile) = data.dgl.dglerr_iris_sfd;
    pltr_my100m_planck(ifile) = data.dgl.onehundomean_planck;
    pltr_my100merr_planck(ifile) = data.dgl.ohmstd_planck;
    pltr_my100merr_beta_planck(ifile) = data.dgl.sem_planck_mc_beta;
    pltr_my100merr_temp_planck(ifile) = data.dgl.sem_planck_mc_temp;
    pltr_my100merr_tau_planck(ifile) = data.dgl.sem_planck_mc_tau;
    % Combine temp and beta errors in quad because they're usually equal in magnitude
    pltr_my100merr_combo_planck(ifile) = sqrt(data.dgl.sem_planck_mc_temp^2 + data.dgl.sem_planck_mc_beta^2);
%     pltr_myohmim_planck(ifile) = data.dgl.ohmim_planck;
    pltr_my100m_iris(ifile) = data.dgl.onehundomean_iris;
    pltr_my100merr_iris(ifile) = data.dgl.ohmstd_iris;
%     pltr_myohmim_iris(ifile) = data.dgl.ohmim_iris;
    pltr_my100m_iris_sfd(ifile) = data.dgl.onehundomean_iris_sfd;
    pltr_my100merr_iris_sfd(ifile) = data.dgl.ohmstd_iris_sfd;
    pltr_mydgl_nh(ifile) = data.dgl.ohmmean_nh;
    pltr_mydglerr_nh(ifile) = data.dgl.ohmstd_nh;
        
    if want_errmags == 1
        myerr(ifile) = data.err_mags.stats.correrr;
    else
        myerr(ifile) = data.stats.correrr;
    end
    
    if want_errmags == 1
        mycrr(ifile) = data.err_mags.ref.bias;
    else
        mycrr(ifile) = data.ref.bias;
    end
    
    if strcmp(flag_method,'new') == 1
        if want_errmags == 1
            myghostdiff(ifile) = data.err_mags.ghost.diffusesub;
            myghostdifferrpos(ifile) = data.err_mags.ghost.diffusesuberrpos;
            myghostdifferrneg(ifile) = data.err_mags.ghost.diffusesuberrneg;
        else
            myghostdiff(ifile) = data.ghost.diffusesub;
            myghostdifferrpos(ifile) = data.ghost.diffusesuberrpos;
            myghostdifferrneg(ifile) = data.ghost.diffusesuberrneg;
%             myghostdiffcomp(ifile) = data.ghost.diffusediff;  
            myscattering_tot(ifile) = data.scattered.gaiasum + data.scattered.masanasum;
            myscattering_gaia(ifile) = data.scattered.gaiasum;
            myscattering_gaia_psferrpos(ifile) = data.scattered.gaiasum_psferrmax - data.scattered.gaiasum;
            myscattering_gaia_psferrneg(ifile) = data.scattered.gaiasum - data.scattered.gaiasum_psferrmin;
            myscattering_gaia_fluxerrpos(ifile) = data.scattered.gaiasum_fluxerrmax - data.scattered.gaiasum;
            myscattering_gaia_fluxerrneg(ifile) = data.scattered.gaiasum - data.scattered.gaiasum_fluxerrmin;
            myscattering_masana(ifile) = data.scattered.masanasum;
            myscattering_masana_psferrpos(ifile) = data.scattered.masanasum_psferrmax - data.scattered.masanasum;
            myscattering_masana_psferrneg(ifile) = data.scattered.masanasum - data.scattered.masanasum_psferrmin;
            myscattering_masana_fluxerrpos(ifile) = data.scattered.masanasum_fluxerrmax - data.scattered.masanasum;
            myscattering_masana_fluxerrneg(ifile) = data.scattered.masanasum - data.scattered.masanasum_fluxerrmin;
        end
    end
    
    mygood = myfieldnum(ifile) == goodfiles;
    
    if sum(mygood) > 0
        % If struct field data.header.bad exists, check if file is good or
        % bad - only mark as good if not bad
        if isfield(data.header,'bad')
            % Only use good images
            if data.header.bad == 0
                isgood(ifile) = 1;
                mystring = sprintf('%f, %d, %s, %f, %f, %f, %d, %f, %f, %f, %f',...
                    mydate(ifile),myfieldnum(ifile),...
                    mytarget{ifile},mysun(ifile),myohm(ifile),myelong(ifile),...
                    mycrr(ifile),...
                    mysig(ifile),myisl(ifile),mydgl(ifile),...
                    mysig(ifile)-myisl(ifile)-mydgl(ifile));
                
                disp(mystring)
                
                image = data.image.calimage;
                save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(ifile),ifile),'image');
                
                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(ifile),ifile),'nanimage');
                
                %figure(1); clf
                %imagesc(data.dgl.dglim)
                %pause(1)
            end
            % If struct field does not exist, files have not been marked good
            % or bad - assume all good
        else
            isgood(ifile) = 1;
                mystring = sprintf('%f, %d, %s, %f, %f, %f, %d, %f, %f, %f, %f',...
                    mydate(ifile),myfieldnum(ifile),...
                    mytarget{ifile},mysun(ifile),myohm(ifile),myelong(ifile),...
                    mycrr(ifile),...
                    mysig(ifile),myisl(ifile),mydgl(ifile),...
                    mysig(ifile)-myisl(ifile)-mydgl(ifile));
                
                disp(mystring)
                
                image = data.image.calimage;
                save(sprintf('../scratch/field%d_image%d.mat',myfieldnum(ifile),ifile),'image');
                
                nanimage = image;
                nanimage(data.mask.onemask) = NaN;
                save(sprintf('../scratch/field%d_masked%d.mat',myfieldnum(ifile),ifile),'nanimage');
                
                %figure(1); clf
                %imagesc(data.dgl.dglim)
                %pause(1)
        end
    end
    
    %if myfieldnum(ifile) == 7
    %  dbstop
    %end
    
end
isgood = logical(isgood);

pltr_thissig = mysig(isgood)-myisl(isgood)-myghostdiff(isgood)-myscattering_tot(isgood);

isgood_masked = mymasked(isgood);
isgood_corrmean = mysig(isgood);
isgood_maskmean = mymaskmean(isgood);
isgood_tri = mytri(isgood);
isgood_psfwing = mypsfwing(isgood);
isgood_isl = myisl(isgood);
isgood_dgl = mydgl(isgood);
isgood_trierr = mytrierr(isgood);
isgood_dglerr = mydglerr(isgood);
if strcmp(flag_method,'new') == 1
    isgood_diffghost = myghostdiff(isgood);
    isgood_diffghosterrpos = myghostdifferrpos(isgood);
    isgood_diffghosterrneg = myghostdifferrneg(isgood);
    isgood_scattering_tot = myscattering_tot(isgood);
    isgood_scattering_gaia = myscattering_gaia(isgood);
    isgood_scattering_gaia_psferrpos = myscattering_gaia_psferrpos(isgood);
    isgood_scattering_gaia_psferrneg = myscattering_gaia_psferrneg(isgood);
    isgood_scattering_gaia_fluxerrpos = myscattering_gaia_fluxerrpos(isgood);
    isgood_scattering_gaia_fluxerrneg = myscattering_gaia_fluxerrneg(isgood);
    isgood_scattering_masana = myscattering_masana(isgood);
    isgood_scattering_masana_psferrpos = myscattering_masana_psferrpos(isgood);
    isgood_scattering_masana_psferrneg = myscattering_masana_psferrneg(isgood);
    isgood_scattering_masana_fluxerrpos = myscattering_masana_fluxerrpos(isgood);
    isgood_scattering_masana_fluxerrneg = myscattering_masana_fluxerrneg(isgood);
end
isgood_field = myfieldnum(isgood);
isgood_b = myb(isgood);

% Calculate mean masked contribution, trilegal contribution, and psf wing
% contribution per field
for i = 1:length(goodfiles)
    k_field = isgood_field == goodfiles(i);
    field_masked = mean(isgood_masked(k_field));
    field_maskmean = mean(isgood_maskmean(k_field));
    field_tri = mean(isgood_tri(k_field));
    field_psfwing = mean(isgood_psfwing(k_field));
    field_corrmean = mean(isgood_corrmean(k_field));
    field_isl = mean(isgood_isl(k_field));
    field_dgl = mean(isgood_dgl(k_field));
    field_trierr = mean(isgood_trierr(k_field));
    field_dglerr = mean(isgood_dglerr(k_field));
    if strcmp(flag_method,'new') == 1
        field_diffghost = mean(isgood_diffghost(k_field));
        field_diffghosterrpos = mean(isgood_diffghosterrpos(k_field));
        field_diffghosterrneg = mean(isgood_diffghosterrneg(k_field));
        field_scattering_tot = mean(isgood_scattering_tot(k_field));
        field_scattering_gaia = mean(isgood_scattering_gaia(k_field));
        field_scattering_gaia_psferrpos = mean(isgood_scattering_gaia_psferrpos(k_field));
        field_scattering_gaia_psferrneg = mean(isgood_scattering_gaia_psferrneg(k_field));
        field_scattering_gaia_fluxerrpos = mean(isgood_scattering_gaia_fluxerrpos(k_field));
        field_scattering_gaia_fluxerrneg = mean(isgood_scattering_gaia_fluxerrneg(k_field));
        field_scattering_masana = mean(isgood_scattering_masana(k_field));
        field_scattering_masana_psferrpos = mean(isgood_scattering_masana_psferrpos(k_field));
        field_scattering_masana_psferrneg = mean(isgood_scattering_masana_psferrneg(k_field));
        field_scattering_masana_fluxerrpos = mean(isgood_scattering_masana_fluxerrpos(k_field));
        field_scattering_masana_fluxerrneg = mean(isgood_scattering_masana_fluxerrneg(k_field));
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
        ' | DGL mean = ',num2str(field_dgl),' | DGL err = ',num2str(field_dglerr),...
        ' | Scattering mean = ',num2str(field_scattering_tot),' | Scattering Gaia = ',num2str(field_scattering_gaia),' | Scattering Masana = ',num2str(field_scattering_masana),...
        ' | Scattering Masana PSF err pos = ',num2str(field_scattering_masana_psferrpos),' | Scattering Masana PSF err neg = ',num2str(field_scattering_masana_psferrneg),' | Scattering Masana Flux err pos = ',num2str(field_scattering_masana_fluxerrpos),' | Scattering Masana Flux err neg = ',num2str(field_scattering_masana_fluxerrneg),...
        ' | Scattering Gaia PSF err pos = ',num2str(field_scattering_gaia_psferrpos),' | Scattering Gaia PSF err neg = ',num2str(field_scattering_gaia_psferrneg),' | Scattering Gaia Flux err pos = ',num2str(field_scattering_gaia_fluxerrpos),' | Scattering Gaia Flux err neg = ',num2str(field_scattering_gaia_fluxerrneg),...
        ' | Diff Ghost mean = ',num2str(field_diffghost),' | Diff Ghost err pos = ',num2str(field_diffghosterrpos),' | Diff Ghost err neg = ',num2str(field_diffghosterrneg)])
end

figure(1); clf
plot(mysun(isgood),mysig(isgood)-myisl(isgood)-mydgl(isgood),'o')
xlabel('Solar Distance')
ylabel('EBL')

if strcmp(flag_method,'new') == 1
    figure(2); clf
    plot(mysun(isgood),myghostdiff(isgood),'o')
    hold on;
    errorbar(mysun(isgood),myghostdiff(isgood), myghostdifferrneg(isgood), myghostdifferrpos(isgood), 'LineStyle','none','Color','k');
    xlabel('Solar Distance')
    ylabel('Summed Ghost Intensity [nW m^{-2} sr^{-1}]')
    
%     figure(3); clf
%     plot(mysun(isgood),myghostdiffcomp(isgood),'o')
%     yline(mean(myghostdiffcomp(isgood)))
%     xlabel('Solar Distance')
%     ylabel('(Predicted Summed Ghost Intensity - ISL Estimation) [nW m^{-2} sr^{-1}]')
end

mydist = zeros(numel(goodfiles),1);
mygal = zeros(numel(goodfiles),1);
mysubmen = zeros(numel(goodfiles),1);
mysuberr = zeros(numel(goodfiles),1);
mymean = zeros(numel(goodfiles),1);
myunc = zeros(numel(goodfiles),1);
mysem = zeros(numel(goodfiles),1);
myavp = zeros(numel(goodfiles),1);
myohmp = zeros(numel(goodfiles),1);
myb_goodfiles = zeros(numel(goodfiles),1);
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

for ifield=1:numel(goodfiles)
    
    whpl = myfieldnum == goodfiles(ifield);
    mydist(ifield) = mean(mysun(whpl & isgood));
    mygal(ifield) = mean(myb(whpl & isgood));
    if strcmp(flag_method,'new') == 1
        thissig = mysig(whpl & isgood)-myisl(whpl & isgood)-mydgl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        pltr_totsig_goodfiles = mysig(whpl & isgood)-myisl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood);
        pltr_thissig_goodfiles(ifield) = mean(mysig(whpl & isgood)-myisl(whpl & isgood)-myghostdiff(whpl & isgood)-myscattering_tot(whpl & isgood));
    else
        thissig = mysig(whpl & isgood)-myisl(whpl & isgood)-mydgl(whpl & isgood);
    end
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
    
    thiserr = myerr(whpl & isgood);
    mysubmen(ifield) = sum(mysig(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    mysuberr(ifield) = std(mysig(whpl & isgood));
    mymean(ifield) = sum(thissig./thiserr.^2)./sum(1./thiserr.^2);
    myunc(ifield) = std(thissig);
    mysem(ifield) = std(thissig)/sqrt(length(thissig));
    pltr_sem(ifield) = std(pltr_totsig_goodfiles)/sqrt(length(pltr_totsig_goodfiles));
    myavp(ifield) = sum(myav(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    myohmp(ifield) = sum(myohm(whpl & isgood)./thiserr.^2)./sum(1./thiserr.^2);
    
end

figure(4); clf
plot(mydist,mysubmen,'o')
hold on
errorbar(mydist,mysubmen,mysuberr,'o');
xlabel('Solar Distance')
ylabel('Error-Weighted Image Mean')

mysubmen

figure(5); clf
%plot(mydist,mymean,'.','MarkerSize',0);
%hold on
errorbar(mydist,mymean,mysem,'.','MarkerSize',20,'MarkerEdge',[0.8500, 0.3250, 0.0980],'LineStyle','none','Color',[0.8500, 0.3250, 0.0980]); %old error bars were myunc, based on std not sem
hold on;
xlabel('Solar Distance')
ylabel('Error-Weighted EBL')

supermean = sum(mymean./mysem.^2)./sum(1./mysem.^2) %old way was myunc instead of mysem
superunc = 1./sqrt(sum(1./mysem.^2)) %old way was myunc instead of mysem

superav = sum(myavp./myunc.^2)./sum(1./myunc.^2)

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

% IRIS error
iris_err = (0.06)/sqrt((1.13*(4.3)^2)/((17.4^2))); % rms noise (0.06 MJy/sr) converted from per iris beam to per lorri image
iris_err_fields = ones(numfields,1)*iris_err;
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
[b_planck, m_planck, sigm_planck, sigb_planck, chi2, q] = fitexy(pltr_my100m_planck_goodfiles.*dl, pltr_thissig_goodfiles, planck_prop_err, pltr_sem);
% iris 
[b_iris, m_iris, sigm_iris, sigb_iris, chi2, q] = fitexy(pltr_my100m_iris_goodfiles.*dl, pltr_thissig_goodfiles, iris_prop_err, pltr_sem);
% iris/sfd 
[b_iris_sfd, m_iris_sfd, sigm_iris_sfd, sigb_iris_sfd, chi2, q] = fitexy(pltr_my100m_iris_sfd_goodfiles.*dl, pltr_thissig_goodfiles, iris_sfd_prop_err, pltr_sem);

p1 = scatter(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,'b');
hold on
p1errx = errorbar(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,pltr_my100merr_combo_planck_goodfiles.*dl,'horizontal', 'LineStyle','none','Color', 'b');
p1erry = errorbar(pltr_my100m_planck_goodfiles.*dl,pltr_thissig_goodfiles,pltr_sem,'vertical', 'LineStyle','none','Color', 'b');
i1 = scatter(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,'g');
i1errx = errorbar(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,iris_err_fields.*dl,'horizontal', 'LineStyle','none','Color', 'g');
i1erry = errorbar(pltr_my100m_iris_goodfiles.*dl,pltr_thissig_goodfiles,pltr_sem,'vertical', 'LineStyle','none','Color', 'g');
i2 = scatter(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,'m');
i2errx = errorbar(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,iris_err_fields.*dl,'horizontal', 'LineStyle','none','Color', 'm');
i2erry = errorbar(pltr_my100m_iris_sfd_goodfiles.*dl,pltr_thissig_goodfiles,pltr_sem,'vertical', 'LineStyle','none','Color', 'm');

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
% hi4pi_err = ((2.3e18)/5)/sqrt((1.13*(16.2)^2)/((17.4^2))); % 5-sig rms sensitivity (2.3e18 cm^-2) converted to 1-sig and from per hi4pi beam to per lorri image
% hi4pi_err_fields = ones(numfields,1)*hi4pi_err;
% n1 = scatter(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,'k');
% hold on
% n1errx = errorbar(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,hi4pi_err_fields,'horizontal', 'LineStyle','none','Color', 'k');
% n1erry = errorbar(pltr_mydgl_nh_goodfiles,pltr_thissig_goodfiles,pltr_sem,'vertical', 'LineStyle','none','Color', 'k');
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

end
