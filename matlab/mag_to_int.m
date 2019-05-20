function int_cnt = mag_to_int( mag, filter )
% Convert from maginute to intensity for values from USNO B1.0

% PARAMETERS
% mag = the magnitute
% 
% int = the intestity

% author: Poppy Immel
% email: pgi8114@rit.edu
% date: 3/14/15 

zero_point = 1;  % need to update this value.

if strcmp( filter, 'B1' ) == 1 | strcmp( filter, 'B2') == 1
	zero_point = 0.438; %lambda_eff for B
else
	zero_point = 0.641; %lambda_eff for R
end

int_cat = zero_point/(10^(0.4*mag));
int_cnt = int_cat;


end

