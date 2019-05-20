%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_photom2cal.m
%%  Jun 2016, MZ
%%  Given an aperture sum in DN S and a magnitude known for that source,
%%  this program computes the conversion factors:
%%  V = -2.5 * log_10(multcal*S) formula
%%  V = 2.5 log_10(S) + lincal formula
%%  and returns both to the user.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [multcal,lincal] = nh_photom2cal(S,mag)

  % this is the multiplicative factor that belongs in the 
  % V = -2.5 * log_10(multcal*S) formula
  multcal = 10.^(-mag/2.5) ./ S;
  
  % this is the additive factor that belongs in the
  % V = 2.5 log_10(S) + V_0 formula
  lincal = -2.5.*log10(multcal);

% end
