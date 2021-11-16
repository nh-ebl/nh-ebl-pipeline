function data = nh_calcdgl(data, paths, flag_method)

dglparams = nh_get_dgl_params();

if strcmp(flag_method,'new') == 1
    %check if planck maps are already made, if not call python script that makes
    %them (isfile returns 1 if file exists)
    if( ~(isfile(sprintf('%splanck_%s_fx.fits',paths.planckdir,data.header.timestamp)) && ...
            isfile(sprintf('%splanck_%s_ra.fits',paths.planckdir,data.header.timestamp)) && ...
            isfile(sprintf('%splanck_%s_dc.fits',paths.planckdir,data.header.timestamp))) )
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
        movefile([pydir,'planck_',data.header.timestamp,'_fx.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_fx.fits']); %move from pydir to paths.planckdir
        movefile([pydir,'planck_',data.header.timestamp,'_ra.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_ra.fits']); %move from pydir to paths.planckdir
        movefile([pydir,'planck_',data.header.timestamp,'_dc.fits'], [paths.planckdir,'planck_',data.header.timestamp,'_dc.fits']); %move from pydir to paths.planckdir
    end
    
    %load in planck maps
    irismap = fitsread(sprintf('%splanck_%s_fx.fits',...
        paths.planckdir,data.header.timestamp));
    irisra = fitsread(sprintf('%splanck_%s_ra.fits',...
        paths.planckdir,data.header.timestamp));
    irisdc = fitsread(sprintf('%splanck_%s_dc.fits',...
        paths.planckdir,data.header.timestamp));
    irisim = irismap;
    
    type = 'planck';
    
elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'old') == 1)
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
end

%Plot DGL image
h = figure(1);
clf;
set(h,'visible','off');
x = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
y = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
imagesc(x,y,(1-((data.dgl.ohmim_iris-dglparams.cib(1))./(data.dgl.ohmim_planck))));
caxis([0 0.4])
axis('xy')
axis image
pbaspect([1 1 1]);
set(gca,'YDir','normal');
xlabel('[arcmin]')
ylabel('[arcmin]')
a = colorbar;
% a.Label.String = '100 \mum Intensity [MJy sr^{-1}]';
a.Label.String = '1 - IRIS/Planck';
% title('Planck');
ext = '.png';
imagename = sprintf('%s%s%s',paths.dgldir,data.header.timestamp,ext);
print(h,imagename, '-dpng');

%calculate mean and std of iris map
ohm_mean = nanmean(irisim(:));
ohm_std = nanstd(irisim(:));

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

if strcmp(type,'iris') == 1
    %nu I_nu is (iris map - 0.8 MJy/sr) * 1e-20 * 3e8/100e-6 * 1e9
    nuinu = (irisim-dglparams.cib(1)) .* dglparams.norm;
elseif strcmp(type,'planck') == 1
    %nu I_nu is (iris map) * 1e-20 * 3e8/100e-6 * 1e9 - planck map already
    %has CIB subtracted
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
dglerr_nuinu = (cbar .* dl).^2 .* dglparams.cib(2).^2;
dglerr_cbar = (nanmean(nuinu(:)) .* dl).^2 .* dglparams.cbar(2).^2;
dglerr_dl = (nanmean(nuinu(:)) .* cbar .* dglparams.A .* 1.1 .* ...
    sqrt(sin(abs(data.coords.galactic(2)).*pi./180))).^2 .* dglparams.g(2).^2;

dgl_err = sqrt(dglerr_nuinu + dglerr_cbar + dglerr_dl);

%whpl = irisim > 10;
%load('lookup/nh_dglcorr.mat');
%corrfac = interp1(dglcorr.ohm,dglcorr.corr,irisim(whpl));
%dglim(whpl) = corrfac .* dglim(whpl);

%mean and std of dgl estimate (at each pixel?)
dgl_mean = nanmean(dglim(:));
dgl_std = nanstd(dglim(:));

if strcmp(type,'iris') == 1
    data.dgl.dglmean_iris = dgl_mean;
    data.dgl.dglerr_iris = dgl_err;
    data.dgl.dglstd_iris = dgl_std;
    data.dgl.ohmmean_iris = ohm_mean;
    data.dgl.ohmstd_iris = ohm_std;
    data.dgl.ohmim_iris = irisim;
    data.dgl.dglim_iris = dglim;
    data.dgl.conv_iris = dglparams.cbar(1).*dl;
    data.dgl.convp_iris = dglparams.norm .* dglparams.cbar(1).*dl;
elseif strcmp(type,'planck') == 1
    data.dgl.dglmean_planck = dgl_mean;
    data.dgl.dglerr_planck = dgl_err;
    data.dgl.dglstd_planck = dgl_std;
    data.dgl.ohmmean_planck = ohm_mean;
    data.dgl.ohmstd_planck = ohm_std;
    data.dgl.ohmim_planck = irisim;
    data.dgl.dglim_planck = dglim;
    data.dgl.conv_planck = dglparams.cbar(1).*dl;
    data.dgl.convp_planck = dglparams.norm .* dglparams.cbar(1).*dl;
end

%Plot DGL image
% h = figure(1);
% clf;
% set(h,'visible','off');
% x = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% y = linspace((-256/2)*4.1/60,(256/2)*4.1/60,257);
% imagesc(x,y,dglim)
% axis('xy')
% axis image
% pbaspect([1 1 1]);
% set(gca,'YDir','normal');
% xlabel('[arcmin]')
% ylabel('[arcmin]')
% a = colorbar;
% a.Label.String = 'DGL Intensity [nW m^{-2} sr^{-1}]';
% title(sprintf('DGL Mean: %.2f',dgl_mean));
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.dgldir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');


%disp(sprintf('Field number: %d; b: %4.2f; DGL: %7.3f; DGLerr: %7.3f.',...
%    data.header.fieldnum,data.coords.galactic(2),dgl_mean,dgl_err))

% end