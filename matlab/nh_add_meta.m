%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_pipeline.m
%%  Jun 2016, MZ
%%  For each file in the NH data directories, this program:
%%   1) Reads it in.
%%   2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data=nh_add_meta(data,flag_method)

A_V = [0.054,0.056,0.220,0.054,0.055,0.155,0.100,3.706,0.811,1.953,...
    2.437,0.305,2.592,0.768];
avfac = 10.^(-A_V./2.5);
A_R = [0.042,0.044,0.174,0.043,0.044,0.123,0.079,2.932,0.6412,1.545,...
    1.928,0.241,2.051,0.608];
arfac = 10.^(-A_R./2.5);

au2km = 1.496e8;

% search for values using keywords
date=data.astrometry.id_spcutcjd;
date_split = strsplit(date,' ');
data.header.date_jd = str2double(date_split(2));
data.header.date_cal = data.astrometry.id_spcutcal;

% extract target names
data.header.target_name = data.astrometry.id_spctcb;

data.header.lambdamin = 350e-9;
data.header.lambdamax = 850e-9;
data.header.lambda = data.astrometry.id_pivot ./ 10;
data.header.launch_jd = 2453755.291667;
data.header.launch_date = '2006-01-19T19:00:00';
data.header.cover_jd = 2453976.500000;
data.header.cover_date = '2006-08-29T00:00:00';
data.header.ccdtemp = data.astrometry.id_ccdtemp;
data.header.exptime = data.astrometry.id_exptime;
% data.header.A_V = A_V(data.header.fieldnum);
% data.header.avfac = avfac(data.header.fieldnum);
% data.header.A_R = A_R(data.header.fieldnum);
% data.header.arfac = arfac(data.header.fieldnum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.cal.nommagoff = data.astrometry.id_photzpt;
if strcmp(flag_method,'new') == 1
    data.cal.magoff = 18.71; %18.71 is lorri-band converted from 18.88 (v-band) from 2020 calibration paper, old 18.53 from fig S5 of nature paper
elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'old') == 1)
    data.cal.magoff = 18.53;
end
data.cal.gaia2lorrimag = -0.0323;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.const.au2km = au2km;
data.const.c = 3e8;
data.const.vegazero = 3636;

% extract distance from target to NH
data.distance.x_NH_target = data.astrometry.id_spctscx;
data.distance.y_NH_target = data.astrometry.id_spctscy;
data.distance.z_NH_target = data.astrometry.id_spctscz;
data.distance.d_NH_target = sqrt(data.distance.x_NH_target.^2 + ...
    data.distance.y_NH_target.^2 + data.distance.z_NH_target.^2) ./ au2km;

% extract distance from target to sun
data.distance.x_sun_target = data.astrometry.id_spctsox;
data.distance.y_sun_target = data.astrometry.id_spctsoy;
data.distance.z_sun_target = data.astrometry.id_spctsoz;
data.distance.d_sun_target = sqrt(data.distance.x_sun_target.^2 + ...
    data.distance.y_sun_target.^2 + ...
    data.distance.z_sun_target.^2) ./ au2km;

% extract distance from NH to sun
data.distance.x_sun_NH = data.astrometry.id_spcsscx;
data.distance.y_sun_NH = data.astrometry.id_spcsscy;
data.distance.z_sun_NH = data.astrometry.id_spcsscz;
data.distance.d_sun_NH = sqrt(data.distance.x_sun_NH.^2 + ...
    data.distance.y_sun_NH.^2 + data.distance.z_sun_NH.^2) ./ au2km;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

galactic = coco([data.astrom.crval1,data.astrom.crval2],...
    'j2000.0','g','d','d');

ecliptic = coco([data.astrom.crval1,data.astrom.crval2],...
    'j2000.0','e','d','d');

data.coords.celestial = [data.astrom.crval1,data.astrom.crval2];
data.coords.galactic = galactic;
data.coords.ecliptic = ecliptic;
data.coords.sol_elon = data.astrometry.id_sol_elon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp_ref = data.ref.eng(:,1:256);

data.ref.engmean = mean(temp_ref(~data.mask.onemask)); % Used in ref corr calculation, doesn't need to change with error unless changing ref corr
if isfield(data.astrom,'biasmthd')
    data.ref.biasmthd = data.astrom.biasmthd; % This is the listed bias method from the header - can be median or mean
    data.ref.biaslevl = data.astrom.biaslevl; % This is the listed bias level from the header
    % Determine currently used bias method
    if strcmp(data.ref.biasmthd,'mean of dark column data') == 1 % If method is mean, add back recorded bias level due to inability to replicate robust mean method
        data.ref.bias = nh_sigclip(data.ref.line) - data.ref.biaslevl; % Calculate new bias level - to be used in nh_calibrate
    elseif strcmp(data.ref.biasmthd,'median of dark column data') == 1
        data.ref.bias = nh_sigclip(data.ref.line) - data.ref.biaslevl; % Calculate new bias level - to be used in nh_calibrate
    else
        fprintf('Absolute panic: bias method not recognized!')
    end
else
    % If old version of old data that doesn't have bias info, subtract
    % median
    data.ref.bias = nh_sigclip(data.ref.line) - median(data.ref.line);
end