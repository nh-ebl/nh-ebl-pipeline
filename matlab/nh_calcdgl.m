function data = nh_calcdgl(data, paths, flag_method)

dglparams = nh_get_dgl_params();

% want = 'planck';
% want = 'planck_mc';
want = 'iris';
% want = 'iris_sfd';
% want = 'nh';

% if strcmp(flag_method,'new') == 1
% Planck map
if strcmp(want,'planck') == 1
    %check if planck maps are already made, if not call python script that makes
    %them (isfile returns 1 if file exists)
    if( ~(isfile(sprintf('%splanck_%s_fx.fits',paths.planckdir,data.header.timestamp))))
        %if at least one of the files isn't there, call python script to make
        %them all
        disp('No Planck file found, retrieving new Planck file.')
        pydir = '/home/symons/nh_ebl_pipeline/py/Planck_Cirrus_Estimation-master/'; %where the python script is
        pyfile = 'get_planck.py'; %name of the python file to run
        imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
        %write the imagefile path in a text file to where the python script is
        fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
        fprintf(fileID,'%s',imagefile); %write imagefile path
        fclose(fileID); %close file
        %call python script
        system(['python ',pydir,pyfile]);
        %get the planck files and move them to their data directory
        if not(isfolder([paths.planckdir]))
            mkdir([paths.planckdir])
        end
        movefile([pydir,'planck_',data.header.timestamp,'_fx.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_fx.fits']); %move from pydir to paths.planckdir
%         movefile([pydir,'planck_',data.header.timestamp,'_ra.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_ra.fits']); %move from pydir to paths.planckdir
%         movefile([pydir,'planck_',data.header.timestamp,'_dc.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_dc.fits']); %move from pydir to paths.planckdir
    end
    
    %load in planck maps
    irismap = fitsread(sprintf('%splanck_%s_fx.fits',...
        paths.planckdir,data.header.timestamp));
%     irisra = fitsread(sprintf('%splanck_%s_ra.fits',...
%         paths.planckdir,data.header.timestamp));
%     irisdc = fitsread(sprintf('%splanck_%s_dc.fits',...
%         paths.planckdir,data.header.timestamp));
    irisim = irismap;
    
    type = 'planck';

    % Planck map
elseif strcmp(want,'planck_mc') == 1
    %check if planck maps are already made, if not call python script that makes
    %them (isfile returns 1 if file exists)
    if( ~(isfile(sprintf('%splanck_%s_errmean.mat',paths.planckmcdir,data.header.timestamp))))
        %if at least one of the files isn't there, call python script to make
        %them all
        disp('No Planck MC file found, retrieving new Planck MC file.')
        pydir = '/home/symons/nh_ebl_pipeline/py/Planck_Cirrus_Estimation-master/'; %where the python script is
        pyfile = 'get_planck_mc.py'; %name of the python file to run
        imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
        %write the imagefile path in a text file to where the python script is
        fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
        fprintf(fileID,'%s',imagefile); %write imagefile path
        fclose(fileID); %close file
        %call python script
        system(['python ',pydir,pyfile]);
        if not(isfolder([paths.planckmcdir]))
            mkdir([paths.planckmcdir])
        end
        %get the planck files and move them to their data directory
        movefile([pydir,'planck_',data.header.timestamp,'_errmean.mat'], [paths.planckmcdir,'planck_',data.header.timestamp,'_errmean.mat']); %move from pydir to paths.planckmcdir
%         movefile([pydir,'planck_',data.header.timestamp,'_ra.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_ra.fits']); %move from pydir to paths.planckdir
%         movefile([pydir,'planck_',data.header.timestamp,'_dc.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_dc.fits']); %move from pydir to paths.planckdir
    end
    
    %load in planck err means
    err_means = load(sprintf('%splanck_%s_errmean.mat',paths.planckmcdir,data.header.timestamp));
    %Tau = 1, Temperature = 2, Beta = 3 | Opacity = 1, Temperature = 2, Spectral-Index = 3
    planck_mc_100m_tau = err_means.err_means(1,:);
    planck_mc_mean_tau = mean(err_means.err_means(1,:));
    planck_mc_sem_tau = std(err_means.err_means(1,:))/sqrt(length(planck_mc_100m_tau));
    planck_mc_100m_temp = err_means.err_means(2,:);
    planck_mc_mean_temp = mean(err_means.err_means(2,:));
    planck_mc_sem_temp = std(err_means.err_means(2,:))/sqrt(length(planck_mc_100m_temp));
    planck_mc_100m_beta = err_means.err_means(3,:);
    planck_mc_mean_beta = mean(err_means.err_means(3,:));
    planck_mc_sem_beta = std(err_means.err_means(3,:))/sqrt(length(planck_mc_100m_beta));

    % Plot histfit and calculated fit together - tau
%     h = figure;
%     clf;
%     set(h,'visible','off');
%     % histfit(maskim(:),bins,'normal');
%     g = histfitCustom(planck_mc_100m_tau,15,'normal');
%     ylims =  ylim;
%     hold on
%     plot(ones(5,1)*planck_mc_mean_tau,linspace(ylims(1),ylims(2),5),'linewidth',2,'color','green');
%     ylim(ylims); % make sure old ylims stay
%     hold off
%     xlabel('100 \mum Image Mean w/ Error in \tau');
%     ylabel('N');
%     title(sprintf('\\mu +/- sem: %.2f +/- %.2f',planck_mc_mean_tau,planck_mc_sem_tau));
%     ext = '.png';
%     type = '_tau';
%     imagename = sprintf('%s%s%s%s',paths.planckmchistdir,data.header.timestamp,type,ext);
%     print(h,imagename, '-dpng');
    
    % Plot histfit and calculated fit together - temp
%     h = figure;
%     clf;
%     set(h,'visible','off');
%     % histfit(maskim(:),bins,'normal');
%     g = histfitCustom(planck_mc_100m_temp,15,'normal');
%     ylims =  ylim;
%     hold on
%     plot(ones(5,1)*planck_mc_mean_temp,linspace(ylims(1),ylims(2),5),'linewidth',2,'color','green');
%     ylim(ylims); % make sure old ylims stay
%     hold off
%     xlabel('100 \mum Image Mean w/ Error in Temp');
%     ylabel('N');
%     title(sprintf('\\mu +/- sem: %.2f +/- %.2f',planck_mc_mean_temp,planck_mc_sem_temp));
%     ext = '.png';
%     type = '_temp';
%     imagename = sprintf('%s%s%s%s',paths.planckmchistdir,data.header.timestamp,type,ext);
%     print(h,imagename, '-dpng');
    
    % Plot histfit and calculated fit together - beta
%     h = figure;
%     clf;
%     set(h,'visible','off');
%     % histfit(maskim(:),bins,'normal');
%     g = histfitCustom(planck_mc_100m_beta,15,'normal');
%     ylims =  ylim;
%     hold on
%     plot(ones(5,1)*planck_mc_mean_beta,linspace(ylims(1),ylims(2),5),'linewidth',2,'color','green');
%     ylim(ylims); % make sure old ylims stay
%     hold off
%     xlabel('100 \mum Image Mean w/ Error in \beta');
%     ylabel('N');
%     title(sprintf('\\mu +/- sem: %.2f +/- %.2f',planck_mc_mean_beta,planck_mc_sem_beta));
%     ext = '.png';
%     type = '_beta';
%     imagename = sprintf('%s%s%s%s',paths.planckmchistdir,data.header.timestamp,type,ext);
%     print(h,imagename, '-dpng');
    
    type = 'planck_mc';

% IRIS map
% elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'old') == 1)
elseif (strcmp(want, 'iris') == 1)
    %       %load in iris maps - old one map per field way
    %       irismap = fitsread(sprintf('%siris_%02d_fx.fits',...
    %           paths.irisdir,data.header.fieldnum));
    %       irisra = fitsread(sprintf('%siris_%02d_ra.fits',...
    %           paths.irisdir,data.header.fieldnum));
    %       irisdc = fitsread(sprintf('%siris_%02d_dc.fits',...
    %           paths.irisdir,data.header.fieldnum));
    %
    %       %create interpolation of iris map and lorri field
    %       F = TriScatteredInterp(irisra(:),irisdc(:),irismap(:));
    %       irisim = F(data.astrometry.ra,data.astrometry.dec);
    
    %check if iris maps are already made, if not call python script that makes
    %them (isfile returns 1 if file exists)
    if( ~(isfile(sprintf('%siris_%s_fx.fits',paths.irisdir,data.header.timestamp)) ))
        %if at least one of the files isn't there, call python script to make
        %them all
        disp('No IRIS file found, retrieving new IRIS file.')
        pydir = '/home/symons/nh_ebl_pipeline/py/IRISpy-master/'; %where the python script is
        pyfile = 'nh_iris_map.py'; %name of the python file to run
        imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
        %write the imagefile path in a text file to where the python script is
        fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
        fprintf(fileID,'%s',imagefile); %write imagefile path
        fclose(fileID); %close file
        %call python script
        system(['python ',pydir,pyfile]);
        if not(isfolder([paths.irisdir]))
            mkdir([paths.irisdir])
        end
        %get the iris files and move them to their data directory
        movefile([pydir,'iris_',data.header.timestamp,'_fx.fits'], [paths.irisdir,'iris_',data.header.timestamp,'_fx.fits']); %move from pydir to paths.irisdir
    end
    %load in iris maps
    irismap = fitsread(sprintf('%siris_%s_fx.fits',...
        paths.irisdir,data.header.timestamp));
%     irisra = fitsread(sprintf('%siris_%s_ra.fits',...
%         paths.irisdir,data.header.timestamp));
%     irisdc = fitsread(sprintf('%siris_%s_dc.fits',...
%         paths.irisdir,data.header.timestamp));
    irisim = irismap;
    
    type = 'iris';

% Combined IRIS/SFD no point-source map
elseif strcmp(want,'iris_sfd') == 1
    %check if iris sfd maps are already made, if not call python script that makes
    %them (isfile returns 1 if file exists)
    if( ~(isfile(sprintf('%siris_sfd_%s_fx.fits',paths.irissfddir,data.header.timestamp)) ))
        %if at least one of the files isn't there, call python script to make
        %them all
        disp('No IRIS SFD file found, retrieving new IRIS SFD file.')
        pydir = '/home/symons/nh_ebl_pipeline/py/IRIS_SFD_Estimation/'; %where the python script is
        pyfile = 'get_iris_sfd.py'; %name of the python file to run
        imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
        %write the imagefile path in a text file to where the python script is
        fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
        fprintf(fileID,'%s',imagefile); %write imagefile path
        fclose(fileID); %close file
        %call python script
        system(['python ',pydir,pyfile]);
        if not(isfolder([paths.irissfddir]))
            mkdir([paths.irissfddir])
        end
        %get the planck files and move them to their data directory
        movefile([pydir,'iris_sfd_',data.header.timestamp,'_fx.fits'], [paths.irissfddir,'iris_sfd_',data.header.timestamp,'_fx.fits']); %move from pydir to paths.irissfddir
%         movefile([pydir,'planck_',data.header.timestamp,'_ra.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_ra.fits']); %move from pydir to paths.planckdir
%         movefile([pydir,'planck_',data.header.timestamp,'_dc.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_dc.fits']); %move from pydir to paths.planckdir
    end
    
    %load in planck maps
    irismap = fitsread(sprintf('%siris_sfd_%s_fx.fits',...
        paths.irissfddir,data.header.timestamp));
%     irisra = fitsread(sprintf('%splanck_%s_ra.fits',...
%         paths.planckdir,data.header.timestamp));
%     irisdc = fitsread(sprintf('%splanck_%s_dc.fits',...
%         paths.planckdir,data.header.timestamp));
    irisim = irismap;
    
    type = 'iris_sfd';

% Neutral Hydrogen map
elseif strcmp(want,'nh') == 1
    %check if planck maps are already made, if not call python script that makes
    %them (isfile returns 1 if file exists)
    if( ~(isfile(sprintf('%snh_%s_fx.fits',paths.nhdir,data.header.timestamp)) ))
        %if at least one of the files isn't there, call python script to make
        %them all
        disp('No NH file found, retrieving new NH file.')
        pydir = '/home/symons/nh_ebl_pipeline/py/NH_Estimation/'; %where the python script is
        pyfile = 'get_nh.py'; %name of the python file to run
        imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
        %write the imagefile path in a text file to where the python script is
        fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
        fprintf(fileID,'%s',imagefile); %write imagefile path
        fclose(fileID); %close file
        %call python script
        system(['python ',pydir,pyfile]);
        if not(isfolder([paths.nhdir]))
            mkdir([paths.nhdir])
        end
        %get the planck files and move them to their data directory
        movefile([pydir,'nh_',data.header.timestamp,'_fx.fits'], [paths.nhdir,'nh_',data.header.timestamp,'_fx.fits']); %move from pydir to paths.nhdir
%         movefile([pydir,'planck_',data.header.timestamp,'_ra.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_ra.fits']); %move from pydir to paths.planckdir
%         movefile([pydir,'planck_',data.header.timestamp,'_dc.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_dc.fits']); %move from pydir to paths.planckdir
    end
    
    %load in planck maps
    irismap = fitsread(sprintf('%snh_%s_fx.fits',...
        paths.nhdir,data.header.timestamp));
%     irisra = fitsread(sprintf('%splanck_%s_ra.fits',...
%         paths.planckdir,data.header.timestamp));
%     irisdc = fitsread(sprintf('%splanck_%s_dc.fits',...
%         paths.planckdir,data.header.timestamp));
    irisim = irismap;
    
    type = 'nh';
end

%Plot DGL image
% h = figure();
% clf;
% set(h,'visible','off');
% x = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% y = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% % imagesc(x,y,(1-((data.dgl.ohmim_iris-dglparams.cib(1))./(data.dgl.ohmim_planck)))); % 1 - IRIS/Planck
% % imagesc(x,y,(irisim-0.48)); % IRIS 100m - CIB
% imagesc(x,y,(irisim)); % 100m or NHI
% % caxis([0 0.4])
% axis('xy')
% axis image
% pbaspect([1 1 1]);
% set(gca,'YDir','normal');
% xlabel('[arcmin]')
% ylabel('[arcmin]')
% a = colorbar;
% a.Label.String = '100 \mum Intensity [MJy sr^{-1}]';
% % a.Label.String = 'Neutral Hydrogen Column Density [cm^{-2}]';
% % a.Label.String = '1 - IRIS/Planck';
% title('IRIS');
% ext = '.png';
% if not(isfolder([paths.dgldir]))
%     mkdir([paths.dgldir])
% end
% imagename = sprintf('%s%s%s',paths.dgldir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

if (strcmp(type,'iris') || strcmp(type,'iris_sfd') || strcmp(type,'nh') || strcmp(type,'planck'))
    if strcmp(type,'iris') == 1
        % Create manual masks
        data = nh_dgl_manmask(data,paths);

        % Load in manual masks
        filein = sprintf('%s%s.mat',paths.dglmandir,data.header.timestamp);
        % If man mask files exists, create the mask, otherwise no man mask
        if numel(dir(filein))
            load(filein);
        else
            manmask = zeros(size(data.data)); % Won't be used - empty masks are saved for all
        end

        %calculate mean and std of iris map
        ohm_mean = mean(mean(irisim(~manmask)));
        ohm_std = std(std(irisim(~manmask)));
    else
        %calculate mean and std of iris map
        ohm_mean = mean(mean(irisim));
        ohm_std = std(std(irisim));
    end

    %[lorri_lambda,lorri_trans] = textread('lookup/nh_lorri_bandpass.txt','%f %f');

    %[dgl_lambda,dgl_conv] = textread('lookup/nh_dgl.txt','%f, %f');

    %dgl_lambda = 1e3 .* dgl_lambda;
    %dgl_conv = dgl_conv ./ 2;

    %lambda = [300:1000];

    %myconv = interp1(dgl_lambda,dgl_conv,lambda,'spline');

    %myconv = sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0) .* ...
    %    myconv) ./ ...
    %    sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0));

    %myfreq = 3e8 ./ 100e-6;

    % If not looking at Neutral Hydrogen, calculate DGL, otherwise skip
    if strcmp(type, 'nh') == 0
        if strcmp(type,'iris') == 1
            %nu I_nu is (iris map - 0.8 MJy/sr) * 1e-20 * 3e8/100e-6 * 1e9
%             nuinu = (irisim.*~manmask-dglparams.cib(1)) .* dglparams.norm;
            onehundo = mean(mean((irisim.*~manmask-0.48))); % mean 100 micron emission including CIB subtraction (real value 0.48, testing 0.24)
            nuinu = (irisim.*~manmask-0.48) .* dglparams.norm; % Old value was 0.8 MJy/sr, now using 0.48 MJy/sr from 14.4 +/- 6.3 nW in https://articles.adsabs.harvard.edu/pdf/2011ASPC..446..309B
        elseif strcmp(type,'iris_sfd') == 1
            %nu I_nu is (iris map - 0.8 MJy/sr) * 1e-20 * 3e8/100e-6 * 1e9
            onehundo = mean(mean(irisim)); % mean 100 micron emission including CIB subtraction
            nuinu = (irisim) .* dglparams.norm; % IRIS/SFD doesn't need CIB sub if zero-referenced to NHI
        elseif strcmp(type,'planck') == 1
            %nu I_nu is (iris map) * 1e-20 * 3e8/100e-6 * 1e9 - planck map already
            %has CIB subtracted
            onehundo = mean(mean(irisim)); % mean 100 micron emission including CIB subtraction
            nuinu = (irisim) .* dglparams.norm;
        end

        %b_lambda (I_nu_opt/I_nu_100um) is estimated by fitting mean ZDA04 model to many measurements
        %of b_lambda
        %c_lambda (lambda*I_lambda_opt/nu*I_nu_100um)s is 10^-6*(100um/0.655um)*b_lambda
        %cbar is 0.491 (cbar_655nm calculated by integrating cbar over lorri
        %bandpass)
        cbar = dglparams.cbar(1);

        %A = 1/0.567 is normalizing factor d0 (from normalizing d(b) at b=25)
        %g = 0.61 is bandpass-weighted mean of observation-constrained model for
        %high-lat diffuse dust component of DGL
        %dl is d(b), a function that accounts for change in c_lambda due to b
        dl = dglparams.A .* (1 - 1.1 .* dglparams.g(1).*sqrt(sin(abs(data.coords.galactic(2)).*pi./180)));

        %dglim is estimate of DGL using average of 100um intensity, c_lambda, and
        %d(b)
        dglim = nuinu .* cbar .* dl;

        %propagated error here
        if strcmp(type,'iris') == 1
            dglerr_nuinu = (cbar .* dl).^2 .* 0.21.^2; % New CIB error is 0.21 MJy/sr
            dglerr_cbar = (mean(nuinu(~manmask)) .* dl).^2 .* dglparams.cbar(2).^2;
            dglerr_dl = (mean(nuinu(~manmask)) .* cbar .* dglparams.A .* 1.1 .* ...
                sqrt(sin(abs(data.coords.galactic(2)).*pi./180))).^2 .* dglparams.g(2).^2;

            dgl_err = sqrt(dglerr_nuinu + dglerr_cbar + dglerr_dl);
        else
            dglerr_nuinu = (cbar .* dl).^2 .* 0.^2; % No error from CIB if not subtracting
            dglerr_cbar = (mean(nuinu(:)) .* dl).^2 .* dglparams.cbar(2).^2;
            dglerr_dl = (mean(nuinu(:)) .* cbar .* dglparams.A .* 1.1 .* ...
                sqrt(sin(abs(data.coords.galactic(2)).*pi./180))).^2 .* dglparams.g(2).^2;

            dgl_err = sqrt(dglerr_nuinu + dglerr_cbar + dglerr_dl);
        end

        %whpl = irisim > 10;
        %load('lookup/nh_dglcorr.mat');
        %corrfac = interp1(dglcorr.ohm,dglcorr.corr,irisim(whpl));
        %dglim(whpl) = corrfac .* dglim(whpl);

        %mean and std of dgl estimate (at each pixel?)
        if strcmp(type,'iris') == 1
            dgl_mean = mean(dglim(~manmask));
            dgl_std = std(dglim(~manmask));
        else
            dgl_mean = mean(dglim(:));
            dgl_std = std(dglim(:));
        end
    end
end

if strcmp(type,'iris') == 1
    data.dgl.dglmean_iris = dgl_mean;
    data.dgl.dglerr_iris = dgl_err;
    data.dgl.dglstd_iris = dgl_std;
    data.dgl.ohmmean_iris = ohm_mean;
    data.dgl.ohmstd_iris = ohm_std;
    data.dgl.ohmim_iris = irisim.*~manmask;
    data.dgl.onehundomean_iris = onehundo;
    data.dgl.dglim_iris = dglim;
    data.dgl.conv_iris = dglparams.cbar(1).*dl;
    data.dgl.convp_iris = dglparams.norm .* dglparams.cbar(1).*dl;
elseif strcmp(type,'iris_sfd') == 1
    data.dgl.dglmean_iris_sfd = dgl_mean;
    data.dgl.dglerr_iris_sfd = dgl_err;
    data.dgl.dglstd_iris_sfd = dgl_std;
    data.dgl.ohmmean_iris_sfd = ohm_mean;
    data.dgl.ohmstd_iris_sfd = ohm_std;
    data.dgl.ohmim_iris_sfd = irisim;
    data.dgl.onehundomean_iris_sfd = onehundo;
    data.dgl.dglim_iris_sfd = dglim;
    data.dgl.conv_iris_sfd = dglparams.cbar(1).*dl;
    data.dgl.convp_iris_sfd = dglparams.norm .* dglparams.cbar(1).*dl;
elseif strcmp(type,'planck') == 1
    data.dgl.dglmean_planck = dgl_mean;
    data.dgl.dglerr_planck = dgl_err;
    data.dgl.dglstd_planck = dgl_std;
    data.dgl.ohmmean_planck = ohm_mean;
    data.dgl.ohmstd_planck = ohm_std;
    data.dgl.ohmim_planck = irisim;
    data.dgl.onehundomean_planck = onehundo;
    data.dgl.dglim_planck = dglim;
    data.dgl.conv_planck = dglparams.cbar(1).*dl;
    data.dgl.convp_planck = dglparams.norm .* dglparams.cbar(1).*dl;
elseif strcmp(type,'nh') == 1
    data.dgl.ohmmean_nh = ohm_mean;
    data.dgl.ohmstd_nh = ohm_std;
    data.dgl.ohmim_nh = irisim;
elseif strcmp(type,'planck_mc') == 1
%     data.dgl.ohm_planck_mc = planck_mc_100m;
%     data.dgl.mean_planck_mc = planck_mc_mean;
%     data.dgl.std_planck_mc = planck_mc_std;
    data.dgl.ohm_planck_mc_tau = planck_mc_100m_tau;
    data.dgl.mean_planck_mc_tau = planck_mc_mean_tau;
    data.dgl.sem_planck_mc_tau = planck_mc_sem_tau;
    data.dgl.ohm_planck_mc_temp = planck_mc_100m_temp;
    data.dgl.mean_planck_mc_temp = planck_mc_mean_temp;
    data.dgl.sem_planck_mc_temp = planck_mc_sem_temp;
    data.dgl.ohm_planck_mc_beta = planck_mc_100m_beta;
    data.dgl.mean_planck_mc_beta = planck_mc_mean_beta;
    data.dgl.sem_planck_mc_beta = planck_mc_sem_beta;
end

% Plot DGL image
% IRIS
% h = figure();
% clf;
% set(h,'visible','off');
% x = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% y = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% imagesc(x,y,data.dgl.dglim_iris)
% axis('xy')
% axis image
% pbaspect([1 1 1]);
% set(gca,'YDir','normal');
% xlabel('[arcmin]')
% ylabel('[arcmin]')
% a = colorbar;
% a.Label.String = 'DGL Intensity [nW m^{-2} sr^{-1}]';
% title(sprintf('IRIS DGL Mean: %.2f',data.dgl.dglmean_iris));
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.dgldir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

% Planck
% h = figure();
% clf;
% set(h,'visible','off');
% x = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% y = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% imagesc(x,y,data.dgl.dglim_planck)
% axis('xy')
% axis image
% pbaspect([1 1 1]);
% set(gca,'YDir','normal');
% xlabel('[arcmin]')
% ylabel('[arcmin]')
% a = colorbar;
% a.Label.String = 'DGL Intensity [nW m^{-2} sr^{-1}]';
% title(sprintf('Planck DGL Mean: %.2f',data.dgl.dglmean_planck));
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.dgldir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');


%disp(sprintf('Field number: %d; b: %4.2f; DGL: %7.3f; DGLerr: %7.3f.',...
%    data.header.fieldnum,data.coords.galactic(2),dgl_mean,dgl_err))

% end