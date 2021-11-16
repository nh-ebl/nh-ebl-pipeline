%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_calibrate.m
%%  Jun 2016, MZ
%%  This function applies known calibration factors to the reducing NH data.
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
data.const.Jy = 1e-26; % W m^-2 Hz^-1 per Jy
data.const.nW = 1e9; % nW per W

% omega of beam compared to one pixel
data.cal.pixperbeam = 2.640; %solid angle of PSF in pix^2

% multiplicative zero point of the conversion
data.cal.szero = 10.^(data.cal.magoff ./ -2.5);

% zero point of the vega system in R_L band
data.cal.vzero = 3050; % Jy at R_L=0.

% compute the number of Jy/bit = 10^-26 W m^-2 Hz^-1 / bit
data.cal.Jyperbit = data.cal.vzero .* data.cal.szero; % Jy/bit

% some astrometric quantities we'll need that are fixed values from elsewhere
if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    data.cal.pixsize = 4.*4.96e-6; % radians, used to be * 1.05 ''/pix -> wrong
elseif strcmp(flag_method, 'old') == 1
    data.cal.pixsize = 4.*4.96e-6*1.05; % radians, used to be * 1.05 ''/pix -> wrong
end
data.cal.pixsize_arcsec = data.cal.pixsize .* 180 .* 3600 ./ pi; % arcsec
data.cal.omega = data.cal.pixperbeam .* data.cal.pixsize.^2;
data.cal.omega_pix = (data.cal.pixsize_arcsec./3600.).^2 * (pi./180.).^2;

% compute the frequency in Hz
data.cal.nu = 3e8 ./ (data.header.lambda .* 1e-9);

% so this becomes the conversion from DN/s to nW m^-2 sr^-1
data.cal.sbconv = data.cal.Jyperbit * data.const.Jy .* data.const.nW .* ...
    data.cal.nu ./ data.cal.omega;

% Old reference correction based on fit from nh_calcref
if strcmp(flag_method, 'old') == 1
    if new == 1
        load('lookup/nh_refcorr_new.mat');
        %     mostprob_date = data.refcorr.mostprob_date;
        %     mostprob_corr = data.refcorr.mostprob_corr;
        mean_date = datastruct.refcorr.date;
        mean_corr = datastruct.refcorr.corr;
        whpl_mean = data.header.date_jd == mean_date;
        %     whpl_mostprob = data.header.date_jd == mostprob_date;
    elseif old == 1
        load('lookup/nh_refcorr_old.mat');
        whpl_mean = data.header.date_jd == correction.date;
        mean_corr = correction.corr;
    elseif ghost == 1
        load('lookup/nh_refcorr_ghost.mat');
        whpl_mean = data.header.date_jd == correction.date;
        mean_corr = correction.corr;
    end
end

% New reference correction
if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    % corr is later subtracted from data, so this becomes med(ref) is added back,
    % sig_mean(ref) is subtracted, and offset is added
    % Offset determined in nh_dark_analysis_combined from robust fit
    newcorr = data.ref.bias - 0.357;
    
    % Adjust data.data directly
    datastruct.data = datastruct.data - newcorr;
    
    % recalculate mask mean and std with new data.data
    datmean = mean(datastruct.data(~datastruct.mask.mask)./ data.astrom.exptime);
    datstd = std(datastruct.data(~datastruct.mask.mask)./ data.astrom.exptime);
    
    % and append it to the data structure
    datastruct.stats.maskmean = datmean;
    datastruct.stats.maskstd = datstd;
    datastruct.stats.maskerr = datstd ./ sqrt(256.^2 - sum(datastruct.mask.onemask(:)));
    
    % calculate some of the stats to take account of the calibration factor
    datastruct.stats.calmean = data.cal.sbconv .* datastruct.stats.maskmean;
    datastruct.stats.calerr = data.cal.sbconv .* datastruct.stats.maskerr;
end

% Reference bias correction used to be added, but seems wrong and now is
% subtracted

% Reference-corrected masked mean in surface brightness units - from saved
% corr file - any old data needs this
% data.stats.corrmean = data.cal.sbconv .* (data.stats.maskmean - ...
%     correction.corr(whpl)./data.header.exptime);
% Reference-corrected masked mean in sb - from data saved corr
if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    % Reference adjustment has already been made directly to data, mean is updated
    datastruct.stats.corrmean = data.cal.sbconv .* (datastruct.stats.maskmean);
elseif strcmp(flag_method, 'old') == 1
    datastruct.stats.corrmean = data.cal.sbconv .* (datastruct.stats.maskmean + ...
        mean_corr(whpl_mean)./data.header.exptime);
end
% Reference-corrected most prob val in sb - from data saved corr
% data.stats.corrmostprob = data.cal.sbconv .* (data.stats.mostprob - ...
%     mostprob_corr(whpl_mostprob)./data.header.exptime);
datastruct.stats.correrr = data.cal.sbconv .* (datastruct.stats.maskerr);

% pull out the raw image data
rawimage = datastruct.data;

% and make calibrated versions of that as well. - if new method, ref corr
% already subtracted from data, if old method, ref corr not applied to calimage
datastruct.image.rawimage = rawimage ./ data.header.exptime;
datastruct.image.calimage = ((rawimage./ data.header.exptime)).* data.cal.sbconv;

if params.err_mags == 1
    data.(params.err_str) = datastruct;
else
    data = datastruct;
end

end


