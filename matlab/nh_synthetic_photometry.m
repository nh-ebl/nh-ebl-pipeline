%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_synthetic_photometry.m
%%  Jun 2016, MZ
%%  This function computes the magnitude of a source in a given optical band
%%  given several USNO-B1 catalog measurements.
%%  Inputs: lambdas = vector of input measurement wavelengths in nm
%%          mags = vector of (Vega) mag measurements at those wavelengths
%%          outband = string defining which optical band to interpolate for
%%  Outputs: mymag = interpolated magnitude at outband
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mymag = nh_synthetic_photometry(lambdas,mags,outband)

  % check we have at least two measurements to work with - this is because
  % if there's only 1 estimate we can't form a combination of it.
  if numel(lambdas) < 2
    warn(sprintf('Too few input data points.'))
  end
  
  % lambda zero points for flux scale
  lambda_zero = [360,438,545,641,798,1220,1630,2190];
  f_zero = [1.79,4.063,3.636,3.064,2.416,1.589,1.021,0.64];
  
  % gridding of lambda out for the magnitude calculation
  fine_lambda = [300:1100];
  
  % interpolate the flux zeros on the new grid
  fine_f_zero = interp1(lambda_zero,f_zero,fine_lambda);
  
  % if there's only two lambdas, make up a third at an intermediate spot so
  % that line fitting will work
  if numel(lambdas) == 2
    lambdas = [lambdas(1),(lambdas(2) + lambdas(1)) ./ 2,lambdas(2)];
    mags = [mags(1),(mags(2) + mags(1)) ./2,mags(2)];
  end
    
  % fit a line to the magnitude measurements to allow interpolation
  fine_mags_coeffs = polyfit(lambdas,mags,1);

  % and use that line model to interpolate
  fine_mags = fine_mags_coeffs(1) .* fine_lambda + fine_mags_coeffs(2);
  
  % if user requested V band fluxes, we need to read in that information
  if strcmp(outband,'V')
    % read in the file
    [filt_lambda,filt_weird,filt_trans] = textread('lookup/vband.txt',...
	'%f %f %f');
    % in the file lambda is quoted as A so convert to nm, and convert from
    % column to row vector
    filt_lambda = filt_lambda' ./ 10;
    % convert from column to row vector
    filt_trans = filt_trans';
  end
  if strcmp(outband,'LORRI')
    % read in the file
    [filt_lambda,filt_trans] = textread('lookup/nh_lorri_bandpass.txt',...
	'%f %f');
    % convert from column to row vector
    filt_lambda = filt_lambda';
    % convert from column to row vector
    filt_trans = filt_trans'./max(filt_trans);
  end
  
  % now interpolate the filter measurements on the same fine grid
  fine_filt_trans = interp1(filt_lambda,filt_trans,fine_lambda,'linear',0.0);
  
  % compute my mag as the filter-weighted integral of the mags
  mymag = sum(fine_mags .* fine_filt_trans) ./ sum(fine_filt_trans);

  % this final catch is to make sure the interpolated magnitude is not
  % wildly different from the input mags, which can happen if we have two
  % very closely spaced lambdas as input
  if abs(mymag - median(mags)) > 2
    mymag = median(mags);
  end

  
% end