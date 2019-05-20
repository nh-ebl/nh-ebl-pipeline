function mylambda = nh_lambda()

  % gridding of lambda out for the magnitude calculation
  fine_lambda = [300:1100];

  % read in the file
  [filt_lambda,filt_trans] = textread('lookup/nh_lorri_bandpass.txt',...
      '%f %f');
  % convert from column to row vector
  filt_lambda = filt_lambda';
  % convert from column to row vector
  filt_trans = filt_trans'./max(filt_trans);

  % now interpolate the filter measurements on the same fine grid
  fine_filt_trans = interp1(filt_lambda,filt_trans,fine_lambda,'linear',0.0);

  mylambda = sum(fine_lambda .* fine_filt_trans) ./ sum(fine_filt_trans);


% end