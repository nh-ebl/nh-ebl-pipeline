%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function wcs_rotate.m
%%  Sept 18, 2010
%%  Mike Zemcov
%%  Rotate between standard (e.g. celestial) and native coordinates
%%  Calling sequence: 
%%       wcsparams = wcs_rotate(wcsparams)
%%  Inputs:
%%       wcsparams = a structure containing the fields lonpole, latpole,
%%                   theta0, crval and either the pair longitude and latitude
%%                   (if you have celestial coordinates) or the pair phi
%%                   and theta (if you have image coordinates). 
%%
%%       Note: If  phi and theta are part of the input parameter structure then
%%       and longitude and latitude are computed.    Otherwise, longitude and
%%       latitude are input structure parameters and phi and theta are computed.
%%
%% Outputs: wcsparams - same as input structure but with fields longitude
%%                      and latitude computed (in the case of theta, phi
%%                      inputs) or phi and theta computed (in the case of
%%                      longitude and latitude as inputs).
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wcsparams=wcs_rotate(wcsparams)

  % define angle constants 
%   pi2 = pi./2.d0;
  radeg = 1.8d2./pi;

  % check that the input fields make sense
  if isfield(wcsparams,'phi') || isfield(wcsparams,'theta')
    if ~isfield(wcsparams,'phi') || ~isfield(wcsparams,'theta')
      fprintf(sprintf(['wcs_rotate error: theta and phi are' ...
		    ' not set correctly!  Aborting.\n']));
      return
    end
    reverseme = 1;
  else 
    if ~isfield(wcsparams,'longitude') || ~isfield(wcsparams,'latitude')
      fprintf(sprintf(['wcs_rotate error: longitude and latitude are' ...
		    ' not set correctly!  Aborting.\n']));
      return
    end
    reverseme = 0;
  end

  if ~isfield(wcsparams,'theta0')
    fprintf(sprintf('wcs_rotate error: theta0 not set!  Aborting.\n'));
  end
  
  % Longpole is the longitude in the native system of the North Pole in the
  % standard system (default = 180 degrees).
  if ~isfield(wcsparams,'lonpole')
    wcsparams.lonpole = 1.8d2;
  end
  phi_p = wcsparams.lonpole ./ radeg;
  sp = sin(phi_p);
  cp = cos(phi_p);

  % If Theta0 = 90 then CRVAL gives the coordinates of the origin in the
  % native system.   This must be converted (using Eq. 7 in Greisen & Calabretta
  % with theta0 = 0) to give the coordinates of the North pole.
  [alpha_p,delta_p] = wcs_getpole(wcsparams);

  % compute useful quantities relating to reference angles
  sa = sin(alpha_p);
  ca = cos(alpha_p);
  sd = sin(delta_p);
  cd = cos(delta_p);

  % calculate rotation matrix 
  r = [ [-sa.*sp - ca.*cp.*sd,  ca.*sp - sa.*cp.*sd, cp.*cd ] ;...
        [ sa.*cp - ca.*sp.*sd, -ca.*cp - sa.*sp.*sd, sp.*cd ] ;...
        [ ca.*cd           ,  sa.*cd           , sd    ] ];
  
  % solve the set of equations for each datum point
  if reverseme
        latitude = wcsparams.phi;
        longitude = wcsparams.theta;
        g = find( isfinite(wcsparams.phi) & isfinite(wcsparams.theta) );
        if numel(g) == 0 
	  % this would be an error
          disp(sprintf('wcs_rotate: No finite input theta, phi given!'));
	  return
	else
	  phi1 = wcsparams.phi(g) ./ radeg;
	  theta1 = wcsparams.theta(g) ./ radeg;
	  r = transpose(r);
	end
  else
    phi = wcsparams.longitude;
    phi1 = wcsparams.longitude ./ radeg;
    theta1 = wcsparams.latitude ./ radeg;
  end

  % define the right-hand side of the equations
  l = cos(theta1) .* cos(phi1);
  m = cos(theta1) .* sin(phi1);
  n = sin(theta1);

 % find solution to the system of equations and put it in b
 % Can't use matrix notation in case l,m,n are vectors
 b0 = r(1,1).*l + r(1,2).*m + r(1,3).*n;
 b1 = r(2,1).*l + r(2,2).*m + r(2,3).*n;
 b2 = r(3,1).*l + r(3,2).*m + r(3,3).*n;

 % use b0,b1,b2 to compute "native" latitude and longitude
 if reverseme
   wcsparams.latitude(g) = asin(b2) .* radeg;
   wcsparams.longitude(g) = atan2( b1, b0) .* radeg;
 else 
   wcsparams.theta = asin(b2) .* radeg;
   wcsparams.phi = atan2( b1, b0) .* radeg;
 end

 return

end
