%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_calibrate.m
%%  Jun 2016, MZ
%%  This function applies known calibration factors to the reducing NH datastruct.
%%  Input = data, an NH pipeline data structrue.
%%  Output = data, an NH pipeline data structure.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = nh_calibrate(data,paths,params,flag_method)

% load('run_params.mat','params')

% set paths here for desired data set
if strcmp(params.data_type,'ghost') == 1
    ghost = 1;
    old = 0;
    new = 0;
elseif strcmp(params.data_type,'old') == 1
    old = 1;
    ghost = 0;
    new = 0;
elseif strcmp(params.data_type,'new') == 1
    new = 1;
    old = 0;
    ghost = 0;
end

if params.err_mags == 1
    datastruct = data.(params.err_str);
else
    datastruct = data;
end

% some constants
datastruct.const.Jy = 1e-26; % W m^-2 Hz^-1 per Jy
datastruct.const.nW = 1e9; % nW per W

% omega of beam compared to one pixel
datastruct.cal.pixperbeam = 2.640; %solid angle of PSF in pix^2

% multiplicative zero point of the conversion
datastruct.cal.szero = 10.^(datastruct.cal.magoff ./ -2.5);

% zero point of the vega system in R_L band
datastruct.cal.vzero = 3050; % Jy at R_L=0.

% compute the number of Jy/bit = 10^-26 W m^-2 Hz^-1 / bit
datastruct.cal.Jyperbit = datastruct.cal.vzero .* datastruct.cal.szero; % Jy/bit

% some astrometric quantities we'll need that are fixed values from elsewhere
if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    datastruct.cal.pixsize = 4.*4.96e-6; % radians, used to be * 1.05 ''/pix -> wrong
elseif strcmp(flag_method, 'old') == 1
    datastruct.cal.pixsize = 4.*4.96e-6*1.05; % radians, used to be * 1.05 ''/pix -> wrong
end
datastruct.cal.pixsize_arcsec = datastruct.cal.pixsize .* 180 .* 3600 ./ pi; % arcsec
datastruct.cal.omega = datastruct.cal.pixperbeam .* datastruct.cal.pixsize.^2;
datastruct.cal.omega_pix = (datastruct.cal.pixsize_arcsec./3600.).^2 * (pi./180.).^2;

% compute the frequency in Hz
datastruct.cal.nu = 3e8 ./ (datastruct.header.lambda .* 1e-9);

% so this becomes the conversion from DN/s to nW m^-2 sr^-1
datastruct.cal.sbconv = datastruct.cal.Jyperbit * datastruct.const.Jy .* datastruct.const.nW .* ...
    datastruct.cal.nu ./ datastruct.cal.omega;

% Old reference correction based on fit from nh_calcref
if strcmp(flag_method, 'old') == 1
    if new == 1
        load('lookup/nh_refcorr_new.mat');
        %     mostprob_date = datastruct.refcorr.mostprob_date;
        %     mostprob_corr = datastruct.refcorr.mostprob_corr;
        mean_date = datastruct.refcorr.date;
        mean_corr = datastruct.refcorr.corr;
        whpl_mean = datastruct.header.date_jd == mean_date;
        %     whpl_mostprob = datastruct.header.date_jd == mostprob_date;
    elseif old == 1
        load('lookup/nh_refcorr_old.mat');
        whpl_mean = datastruct.header.date_jd == correction.date;
        mean_corr = correction.corr;
    elseif ghost == 1
        load('lookup/nh_refcorr_ghost.mat');
        whpl_mean = datastruct.header.date_jd == correction.date;
        mean_corr = correction.corr;
    end
end

% New reference correction
if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    % corr is later subtracted from data, so this becomes med(ref) is added back,
    % sig_mean(ref) is subtracted, and offset is added
        
    % Offset determined in nh_dark_analysis_combined from robust fit
    newcorr = datastruct.ref.bias - 0.357;
    
    % Apply correction to jail bar-removed data and replace data.data
    datastruct.data = datastruct.image.jailbarrem_data - newcorr;
        
    % recalculate mask mean and std with new datastruct.data
    datmean = mean(datastruct.data(~datastruct.mask.mask)./ datastruct.astrom.exptime);
    datstd = std(datastruct.data(~datastruct.mask.mask)./ datastruct.astrom.exptime);
    
    % and append it to the data structure
    datastruct.stats.maskmean = datmean; % This is the mean of the reference-corrected image in DN/s
    datastruct.stats.maskmean_dn = datmean.*datastruct.astrom.exptime; % This is the mean of the reference-corrected image in DN
    datastruct.stats.maskstd = datstd;
    datastruct.stats.maskerr = datstd ./ sqrt(256.^2 - sum(datastruct.mask.onemask(:)));
    
    % calculate some of the stats to take account of the calibration factor
    datastruct.stats.calmean = datastruct.cal.sbconv .* datastruct.stats.maskmean;
    datastruct.stats.calerr = datastruct.cal.sbconv .* datastruct.stats.maskerr;
end

% Reference bias correction used to be added, but seems wrong and now is
% subtracted

% Reference-corrected masked mean in surface brightness units - from saved
% corr file - any old data needs this
% datastruct.stats.corrmean = datastruct.cal.sbconv .* (datastruct.stats.maskmean - ...
%     correction.corr(whpl)./datastruct.header.exptime);
% Reference-corrected masked mean in sb - from data saved corr
if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    % Reference adjustment has already been made directly to data, mean is updated
    datastruct.stats.corrmean = datastruct.cal.sbconv .* (datastruct.stats.maskmean);
elseif strcmp(flag_method, 'old') == 1
    datastruct.stats.corrmean = datastruct.cal.sbconv .* (datastruct.stats.maskmean + ...
        mean_corr(whpl_mean)./datastruct.header.exptime);
end
% Reference-corrected most prob val in sb - from data saved corr
% datastruct.stats.corrmostprob = datastruct.cal.sbconv .* (datastruct.stats.mostprob - ...
%     mostprob_corr(whpl_mostprob)./datastruct.header.exptime);
datastruct.stats.correrr = datastruct.cal.sbconv .* (datastruct.stats.maskerr);

% pull out the raw image data
rawimage = datastruct.data;

% and make calibrated versions of that as well. - if new method, ref corr
% already subtracted from data, if old method, ref corr not applied to calimage
datastruct.image.rawimage = rawimage ./ datastruct.header.exptime; % This is the reference-corrected data in DN/s
datastruct.image.calimage = ((rawimage./ datastruct.header.exptime)).* datastruct.cal.sbconv; % This is the reference-corrected data in nW

if params.err_mags == 1
    data.(params.err_str) = datastruct;
else
    data = datastruct;
end

%Plot masked data with optional mouse movement returns value
% h = figure(1);
% clf;
% imagesc(datastruct.data.*~datastruct.mask.onemask)
% set(h,'visible','off');
% % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
% a = colorbar;
% a.Label.String = 'Intensity [DN]';
% pbaspect([1 1 1]);
% xlabel('LORRI X Pixels');
% ylabel('LORRI Y Pixels');
% caxis([-10,10]);
% %         title(sprintf('%s',datastruct.header.rawfile));
% % grid minor;
% % title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
% title(sprintf('Field: %d',datastruct.header.fieldnum));
% set(gca,'YDir','normal');
% ext = '.png';
% if not(isfolder([paths.corrmaskdir]))
%     mkdir([paths.corrmaskdir])
% end
% imagename = sprintf('%s%s%s',paths.corrmaskdir,datastruct.header.timestamp,ext);
% % imagename = sprintf('%s%s%s',paths.selectmaskdir,datastruct.header.timestamp,ext);
% print(h,imagename, '-dpng');

end


