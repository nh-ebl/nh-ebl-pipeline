function [alpha_p,delta_p] = wcs_getpole(wcsparams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function wcs_getpole.m
%%  Sept. 21, 2010
%%  Mike Zemcov
%%
%%  This function determines the native pole of the coordinate system for
%%  measuring astrometry.
%%
%%  Calling sequence: 
%%       wcsparams = wcs_rotate(wcsparams)
%%
%%  Inputs:
%%       wcsparams = a structure containing the fields lonpole, latpole,
%%                   theta0, and crval.
%%
%%  Outputs: alpha_p     - ra corresponding to coordinate system pole
%%           delta_p     - dec corresponding to coordinate system pole
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% DEFINE ANGLE CONSTANTS 
  pi2 = pi ./ 2.0d0;
  radeg = 1.8d2 ./ pi;
  alpha_0 = wcsparams.crval(1) ./ radeg;
  delta_0 = wcsparams.crval(2) ./ radeg;

  %% if theta0 is 90 then this is simple
  if wcsparams.theta0 == 90
    alpha_p = alpha_0;
    delta_p = delta_0;
    return
  end
  
  %% Longpole is the longitude in the native system of the North Pole in the
  %% standard system (default = 180 degrees).
  phi_p = wcsparams.lonpole ./ radeg;
  theta_p = wcsparams.latpole ./ radeg;
  sp = sin(phi_p);
  cp = cos(phi_p);
  sd = sin(delta_0);
  cd = cos(delta_0);
  tand = tan(delta_0);

  %% If origin is set then crval gives the coordinates of the origin in the
  %% native system.   This must be converted (using Eq. 7 in 
  %% Greisen & Calabretta with theta0 = 0) to give the coordinates
  %% of the North pole (alpha_p, delta_p) 
  if (wcsparams.theta0 == 0)  %% for theta0=0
    if ( (delta_0 == 0) & (abs(wcsparams.lonpole) == 90.0d0) )
      delta_p = theta_p; 
    else 
      delta_p = acos(sd./cp);
    end
    if ( wcsparams.latpole ~= 90 & ...
        abs(theta_p + delta_p) < abs(theta_p - delta_p) )
	delta_p = -delta_p; 
    end
    if ( (wcsparams.lonpole == 1.8d2) || (cd == 0.0) )
      alpha_p = alpha_0; 
    else 
      alpha_p = alpha_0 - atan2(sp./cd, -tan(delta_p).*tand );
    end
  else                 %%General case for arbitary theta0 - if you're in
                       %here it's not going to get any more fun
    ctheta = cos(wcsparams.theta0 ./ radeg);
    stheta = sin(wcsparams.theta0 ./ radeg);
    term1 = atan2(stheta, ctheta .* cp ); 
    term2 = acos( sd ./ ( sqrt(1.0 - ctheta.^2 .* sp.^2)  ));
    if term2 == 0.0
      delta_p = term1; 
    else 
      delta_p1 = abs( (term1 + term2) .* radeg);
      delta_p2 = abs( (term1 - term2) .* radeg);
	
      if (delta_p1 > 90) & (delta_p2 > 90)
        disp(sprintf('wcs_getpole: No valid solution!'));
      end
      if (delta_p1 <= 90) & (delta_p2 > 90)
        delta_p = term1 + term2;
      end
      if (delta_p1 > 90) & (delta_p2 <= 90)
        delta_p = term1 - term2;
      end
      if (delta_p1 <= 90) & (delta_p2 <= 90)
        %%Oh no, there are two valid solutions
        delta_p1 = (term1 + term2) .* radeg;
        delta_p2 = (term1 - term2) .* radeg;
        if abs(wcsparams.latpole - delta_p1) < ... 
              abs(wcsparams.latpole - delta_p2)
          delta_p = term1 + term2;
        else 
          delta_p = term1 - term2;
        end
      end
	
      if (cd == 0.0) 
        alpha_p = alpha_0;
      else
        sdelt = sin(delta_p);
	if (sdelt == 1) 
          alpha_p = alpha_0 - phi_p - pi;
        else 
          if (sdelt == -1) 
            alpha_p = alpha_0 -phi_p; 
          else 
            alpha_p = alpha_0 - ...
                      atan2( (stheta - sin(delta_p).*sd) ./ ...
                             (cos(delta_p).*cd), sp.*ctheta./cd );
          end
        end
      end
      
    end
  end
  
  return
 
end
