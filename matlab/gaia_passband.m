clear all
close all

% read in the LORRI file
[filt_lambda,filt_trans] = textread('lookup/nh_lorri_bandpass.txt','%f %f');

% read in the Gaia file
[gaia_lambda,gaia_trans,null,null,null,null,null] = textread('lookup/GaiaDR2_Passbands.dat','%f %f %f %f %f %f %f');

% Plot original

figure(1)
plot(filt_lambda,filt_trans)
hold on
plot(gaia_lambda(gaia_trans~=99.99),gaia_trans(gaia_trans~=99.99))
xlabel('Wavelength [nm]')
ylabel('Transmissivity')

% gaia_trans(gaia_trans~=99.99) = gaia_trans(gaia_trans~=99.99)./max(gaia_trans(gaia_trans~=99.99)); % Scale gaia_trans to max of 1

% Plot filt scaled
figure(2)
plot(filt_lambda,filt_trans./max(filt_trans))
hold on
plot(gaia_lambda(gaia_trans~=99.99),gaia_trans(gaia_trans~=99.99))
xlabel('Wavelength [nm]')
ylabel('Transmissivity')

lorri_pivot_trans = interp1(filt_lambda,filt_trans./max(filt_trans),607.62); % Get lorri_pivot_trans value at pivot lambda

gaia_pivot_trans = interp1(gaia_lambda,gaia_trans,607.62); % Get gaia_pivot_trans value at pivot lambda

gaia_trans_scaled = gaia_trans./gaia_pivot_trans.*lorri_pivot_trans; % Scale gaia_trans to lorri_pivot_trans value at pivot lambda

% Plot filt and gaia scaled
figure(3)
plot(filt_lambda,filt_trans./max(filt_trans))
hold on
plot(gaia_lambda(gaia_trans~=99.99),gaia_trans_scaled(gaia_trans~=99.99))
xlabel('Wavelength [nm]')
ylabel('Transmissivity')

filt_integral = trapz(filt_lambda,filt_trans./max(filt_trans)); % Integrate under the pivot-wavelength-scaled curve
gaia_integral = trapz(gaia_lambda(gaia_trans~=99.99),gaia_trans_scaled(gaia_trans~=99.99)); % Integrate under the pivot-wavelength-scaled curve

filt_per_gaia_integral = filt_integral/gaia_integral; % Calculate the ratio between the integrals
gaia_per_filt_integral = gaia_integral/filt_integral; % Calculate the ratio between the integrals

disp(['Filt/gaia integral ratio: ',num2str(filt_per_gaia_integral),' | Gaia/filt integral ratio: ',num2str(gaia_per_filt_integral)])

fprintf('done');