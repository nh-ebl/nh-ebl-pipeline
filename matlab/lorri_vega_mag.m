% Calculate LORRI-band magnitude of Vega from other bands

% Lambdas and magnitudes taken from Table 2 of Megessier 1995 (UBVRIJ)
% Different answers whether using Johnson 1964, 1965, or 1966 values
lambda = [360,440,550,700,900,1250];
mags = [0.03,0.03,0.03,0.07,0.10,0.02];

% Output from synthetic photometry (bandpass-weighted)
vega_Rl = nh_synthetic_photometry(lambda,mags,'LORRI');
display(vega_Rl)