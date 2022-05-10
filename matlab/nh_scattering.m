% This function calculates the contribution from extended scattering of
% sources due to the extended PSF from Lauer

% A python script is called to calculate the contribution from the Masana
% ISL map beyond 5 degrees, if saved files do not already exist

% Sources from Gaia between the ghost radius and 5 degrees are estimated
% separately here, and then the combined sum is taken as the total
% contribution

function data = nh_scattering(data, paths, errflag_mags, params)

% Read in results of python function here
%check if masana files are already made, if not call python script that makes
%them (isfile returns 1 if file exists)
if( ~(isfile(sprintf('%s%s.mat',paths.scatteringdir,data.header.timestamp))))
    %if at least one of the files isn't there, call python script to make
    %them all
    disp('No scattering file found, creating new scattering file.')
    pydir = '/home/symons/nh_ebl_pipeline/py/scattering/'; %where the python script is
    pyfile = 'get_scattering_image.py'; %name of the python file to run
    imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
    %write the imagefile path in a text file to where the python script is
    fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
    fprintf(fileID,'%s',imagefile); %write imagefile path
    fclose(fileID); %close file
    %call python script
    system(['python ',pydir,pyfile]);
    %get the scattering files and move them to their data directory
    if not(isfolder([paths.scatteringdir]))
        mkdir([paths.scatteringdir])
    end
    movefile([pydir,data.header.timestamp,'.mat'], [paths.scatteringdir,data.header.timestamp,'.mat']); %move from pydir to paths.scatteringdir
end

% Read in saved masana info from python script
masana = load(sprintf('%s%s.mat',paths.scatteringdir,data.header.timestamp));
masana_scattered = masana.scattered_nW;
masana_radbins = masana.rad_bins_cent;
masana_psfmaxerr = masana.psf_max_err;
masana_psfminerr = masana.psf_min_err;
masana_fluxmaxerr = masana.flux_max_err;
masana_fluxminerr = masana.flux_min_err;
masana_tot = sum(masana_scattered);
masana_tot_psfmaxerr = sum(masana_psfmaxerr);
masana_tot_psfminerr = sum(masana_psfminerr);
masana_tot_fluxmaxerr = sum(masana_fluxmaxerr);
masana_tot_fluxminerr = sum(masana_fluxminerr);

% load lauer psf
lauer_psf = readmatrix('lookup/lauer_psf.csv');
% ang_dist_big = linspace(0.3037,4,100);
% bright_big = interp1(lauer_psf(:,1),lauer_psf(:,2),ang_dist_big,'spline');

% Plot interp check
% figure(1)
% plot(log10(lauer_psf(:,1),lauer_psf(:,2),ang_dist_big,'spline');
% hold('on');
% plot(log10(ang_dist_big),log10(bright_big));

% load up the corresponding catalog file
load(sprintf('%smat_files/field_%d_data_wide.mat',paths.gaiadir,data.header.fieldnum));
Gmag_par = Gmag;
RA_par = RA;
DEC_par = DEC;

% figure out the length of the catalog
[ncat,~] = size(RA);

% preallocate
raa=zeros(ncat,1);
decc=zeros(ncat,1);
boxxpix = zeros(ncat,1);
boxypix = zeros(ncat,1);
boxmag = zeros(ncat,1);
boxmag_errmax = zeros(ncat,1);
boxmag_errmin = zeros(ncat,1);
stardistcentbox = zeros(ncat,1);

% loop over each catalog entry;
if isempty(gcp('nocreate')) % Check if parpool is active
    maxNumCompThreads(7); % Declare number of threads to use
    parpool('threads'); % Fire up threaded parfor for 7 threads
end
astrom = data.astrom; % Create variable for data.astrom so parfor doesn't copy data to each core's process
parfor row = 1:ncat % Parallel for

    %mags don't seem to have any wild values, so we'll use all of
    %them without restriction
    thismag = Gmag_par(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value

    %find x/y coordinate of the object
    [ypix, xpix] = radec2pix(RA_par(row),DEC_par(row), astrom);

    % Calculate distance from center pixel to star pixel and
    % distance from star pixel to ghost pixel
    stardistcentbox(row) = sqrt((xpix-(128)).^2 + (ypix-(128)).^2);

    % If star between ghost range 18.67 arcmin (274.5 pix) and 5 degrees
    % (4411.76 pix) - 4.08''/pix
    % include in calculation
    if stardistcentbox(row) > 274.5 && stardistcentbox(row) <= 4411.76
        boxxpix(row) = xpix;
        boxypix(row) = ypix;
        boxmag(row) = thismag;
    end
end

% remove zero values
boxxpix = boxxpix(boxxpix~=0);
boxypix = boxypix(boxypix~=0);
stardistcentbox = stardistcentbox(boxmag~=0);
boxmag = boxmag(boxmag~=0);

% Calculate stars individually

%Calculate flux/vega for star
F_box = 10.^(-boxmag/2.5);
F_box_errmax = (10.^(-boxmag/2.5)) + (10.^(-boxmag/2.5)).*0.0017; %Get max error on flux based on converting zero pt mag error to flux
F_box_errmin = (10.^(-boxmag/2.5)) - (10.^(-boxmag/2.5)).*0.0017; %Get min error on flux

% Calculate lauer psf value at star's angle
star_ang = stardistcentbox*4.08/3600; % Convert radius in pixels to degrees
psf_star = interp1(lauer_psf(:,1),lauer_psf(:,2),star_ang,'spline'); % Calculate PSF value at that angle
psf_star_max = psf_star + (psf_star.*0.1); % Maximum with 10% error
psf_star_min = psf_star - (psf_star.*0.1); % Minimum with 10% error

% Scale flux bins by psf beam in DN s^-1 pix^-1
L_G_scattered = F_box.*psf_star;
L_G_scattered_fluxerrmax = F_box_errmax.*psf_star;
L_G_scattered_fluxerrmin = F_box_errmin.*psf_star;
L_G_scattered_psferrmax = F_box.*psf_star_max;
L_G_scattered_psferrmin = F_box.*psf_star_min;
L_G_scattered_nW = L_G_scattered.*406.5386; % Convert from DN/s to nW m^-2 sr^-1
L_G_scattered_nW_fluxerrmax = L_G_scattered_fluxerrmax.*406.5386; % Convert from DN/s to nW m^-2 sr^-1
L_G_scattered_nW_fluxerrmin = L_G_scattered_fluxerrmin.*406.5386; % Convert from DN/s to nW m^-2 sr^-1
L_G_scattered_nW_psferrmax = L_G_scattered_psferrmax.*406.5386; % Convert from DN/s to nW m^-2 sr^-1
L_G_scattered_nW_psferrmin = L_G_scattered_psferrmin.*406.5386; % Convert from DN/s to nW m^-2 sr^-1

% Calculate total summed contribution from gaia
gaia_tot = sum(L_G_scattered_nW);
gaia_tot_psferrmax = sum(L_G_scattered_nW_psferrmax);
gaia_tot_psferrmin = sum(L_G_scattered_nW_psferrmin);
gaia_tot_fluxerrmax = sum(L_G_scattered_nW_fluxerrmax);
gaia_tot_fluxerrmin = sum(L_G_scattered_nW_fluxerrmin);

% Calculate stars in bins
% % preallocate
% num_mag = zeros((length(bright_big)-1),1);
% lower_ang = zeros((length(bright_big)-1),1);
% upper_ang = zeros((length(bright_big)-1),1);
% for i = 1:(length(bright_big)-1)
%     % All star mags that fall into bin
%     mag_bin = boxmag(stardistcentbox >= ang_dist_big(i) && stardistcentbox < ang_dist_big(i+1));
%
%     % Save number of stars in bin and lower and upper angle of bin
%     num_mag(i) = length(mag_bin);
%     lower_ang(i) = ang_dist_big(i);
%     upper_ang(i) = ang_dist_big(i+1);
%
%     % Calculate total lIl in bin
%     bin_area = pi*(upper_ang(i)^2 - lower_ang(i)^2);
%     omega = bin_area*(pi/180)^2;
%     F_bin = data.cal.vzero .* 10.^(-mag_bin/2.5);
%     lIl_bin = 1e-26.*1e9.*data.cal.nu.*F_bin./(omega);
%
%     % Calculate Vega flux in W m^-2 at LORRI wavelength
%     vega_flux = data.cal.vzero*(data.const.Jy)*data.cal.nu; % Vega flux for LORRI band converted from Jy to W m^-2
%
%     % Scale flux bins by psf beam in DN s^-1 pix^-1
%     L_G_scattered = L_G_rad/vega_flux*psf_rad_bins;
%     L_G_scattered_new = L_G_rad_new/vega_flux*psf_rad_bins;
%     % L_G_scattered = L_G_vegas*bright_range
%     % L_G_scattered_new = L_G_vegas_new*bright_range
%     L_G_scattered_nW = L_G_scattered*406.5386; % Convert from DN/s to nW m^-2 sr^-1
%     L_G_scattered_new_nW = L_G_scattered_new*406.5386;
%
%     % Scale lIl by lauer psf
% end

% If error is on, write to error substruct
% load('run_params.mat','params')
if params.err_mags == 1
    % save scattering info to err struct
    data.(params.err_str).scattered.gaiasum = gaia_tot;
    data.(params.err_str).scattered.masanasum = masana_tot;
else
    % save scattering info to data
    data.scattered.gaiasum = gaia_tot;
    data.scattered.gaiasum_psferrmax = gaia_tot_psferrmax;
    data.scattered.gaiasum_psferrmin = gaia_tot_psferrmin;
    data.scattered.gaiasum_fluxerrmax = gaia_tot_fluxerrmax;
    data.scattered.gaiasum_fluxerrmin = gaia_tot_fluxerrmin;
    data.scattered.masanasum = masana_tot;
    data.scattered.masanasum_psferrmax = masana_tot_psfmaxerr;
    data.scattered.masanasum_psferrmin = masana_tot_psfminerr;
    data.scattered.masanasum_fluxerrmax = masana_tot_fluxmaxerr;
    data.scattered.masanasum_fluxerrmin = masana_tot_fluxminerr;
    data.scattered.masana_radbins = masana_radbins;
end
end