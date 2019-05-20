function params = nh_get_dgl_params()

  params.norm = 1e-20 .* 3e8./100e-6 .* 1e9;
  params.cib = [0.80,0.25];
  params.cbar = [0.491,0.129];
  % from Draine; from Sano
  params.g = [0.61,0.10];%[0.752,0.098];
  % from Draine/Lillie at 25 deg; from Sano at 40 degrees
  params.A = 1 ./ 0.567;%(0.35);

% end