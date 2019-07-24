function [xpix,ypix]=wcsxy2sph(longitude,latitude,wcsparams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function wcssph2xy.m
%%  Sept. 21, 2010
%%  Mike Zemcov
%%
%%  This function takes a set ra,dec pairs and a stripped down version of
%%  the astrometry structure and outputs the corresponding pixel values
%%  on a grid.  It follows the instructions in Calabretta & Greisen,
%%  A&A 395, 1077 (2002). 
%%
%%  Inputs: longitude - ra input vector
%%           latitude - dec input vector
%%          wcsparams - a stripped down fits astrometry structure
%%                      containing lonpole, latpole, and crval from the header.
%%
%%  Outputs: x - a list of zeroed and scaled x angles corresponding to
%%               the input ra list
%%           y - a list of zeroed ans scaled y angles corresponding to
%%               the input dec list
%% 
%%  Notes: Right now, this function only supports the TAN projection;
%%         more can be added but this requires more complexity in the code.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %% define angle constants 
   pi2 = pi ./ 2;
   radeg = 180 ./ pi;
   
   %% find the number of elements in each of the data arrays
   n_long = numel( longitude );
   n_lat  = numel( latitude );

   %% error checking
   if (n_long ~= n_lat)
     disp(sprintf(['wcssph2xy.m: longitude and latitude inputs must'...
                   'have the same number of elements.'])); 
   end
   if ~isfield(wcsparams,'lonpole')
     disp(sprintf(['wcssph2xy.m: wcsparams structure field lonpole ' ...
                   'missing!']));
     return
   end
   if ~isfield(wcsparams,'latpole')
     disp(sprintf(['wcssph2xy.m: wcsparams structure field latpole ' ...
                   'missing!']));
     return
   end 
   if ~isfield(wcsparams,'crval')
     disp(sprintf(['wcssph2xy.m: wcsparams structure field crval ' ...
                   'missing!']));
     return
   else
     if (numel(wcsparams.crval) ~= 2)
       disp(sprintf(['wcssph2xy.m: wcsparams structre field crval must ' ...
                     'have 2 elements!']));
       return
     end
   end

   %% Convert all longitude values into the range -180 to 180 so that equations
   %% work properly.
   lng = longitude;
   lat = latitude;
   temp = find(lng > 180.0);
   Ntemp = numel(temp);
   if Ntemp > 0 
     lng(temp) = lng(temp) - 360.;
   end

   %% Make small offsets at poles to allow the transformations to be 
   %% completely invertible.  These generally introduce a small fractional
   %% error  but only at the poles.  They are necessary since all maps  lose
   %% information at the poles when a rotation is applied, because all
   %% points within offset of the poles are mapped to  the same points.
   offset = 1.d-7;
   bad = find(abs(lat - 90.0) < offset*radeg);
   if numel(bad) > 0 
     lat(bad) = 90.0 - north_offset .* radeg;
   end
   bad = find(abs(lat + 90.0) < offset*radeg);
   if numel(bad) > 0 
     lat(bad) = south_offset*radeg - 90.0;
   end

   %% append to wcsparams for passing to rotation function
   wcsparams.longitude = lng;
   wcsparams.latitude  = lat;
   wcsparams.theta0    = 90.0;

   %% Convert from standard coordinate system to "native" coordinate system
   wcsparams = wcs_rotate(wcsparams);

   %% strip out the parts of interest
   phi = wcsparams.phi ./ radeg;
   theta = wcsparams.theta ./ radeg;

   %% This inversion is only for the TAN projection - more can be added
   %% but that requires this get changed.
   xpix = NaN .* zeros(size(theta));
   ypix = xpix;
   g = find(theta > 0);
   Ng = numel(g);
   if Ng > 0
     r_theta = radeg ./ tan(theta(g));
     xpix(g) = r_theta .* sin(phi(g));
     ypix(g) = -r_theta .* cos(phi(g));
   end
   
   %% and that's all she wrote
   return
   
end
