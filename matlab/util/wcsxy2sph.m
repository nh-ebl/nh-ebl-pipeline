function [longitude,latitude]=wcsxy2sph(x,y,wcsparams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function wcsxy2sph.m
%%  Sept. 14, 2010
%%  Mike Zemcov
%%
%%  This function takes a set pixel pairs and a stripped down version of
%%  the astrometry structure and outputs the corresponding ra and dec of
%%  those pixels.  It follows the instructions in Calabretta & Greisen,
%%  A&A 395, 1077 (2002). 
%%
%%  Inputs: x - a list of zeroed and scaled x angles
%%          y - a list of zeroed ans scaled y angles
%%          wcsparams - a stripped down fits astrometry structure
%%                      containing lonpole, latpole, and crval from the header.
%%
%%  Outputs: longitude - ra corresponding to the input x vector
%%           latitude  - dec corresponding to the input y vector
%% 
%%  Notes: Right now, this function only supports the TAN projection;
%%         more can be added but this requires more complexity in the code.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   %% Define angle constants 
   radeg = 180.d0./pi;
   pi2 = pi/2.0d0;

   %% GENERAL ERROR CHECKING
   if ~isfield(wcsparams,'lonpole')
     disp(sprintf(['wcsxy2sph.m: wcsparams structure field lonpole ' ...
                   'missing!']));
     return
   end
   if ~isfield(wcsparams,'latpole')
     disp(sprintf(['wcsxy2sph.m: wcsparams structure field latpole ' ...
                   'missing!']));
     return
   end 
   if ~isfield(wcsparams,'crval')
     disp(sprintf(['wcsxy2sph.m: wcsparams structure field crval ' ...
                   'missing!']));
     return
   else
     if (numel(wcsparams.crval) ~= 2)
       disp(sprintf(['wcsxy2sph.m: wcsparams structre field crval must ' ...
                     'have 2 elements!']));
       return
     end
   end
   
   %% find the number of elements in each of the data arrays
   n_x = numel(x);
   n_y = numel(y);
   sz_x = size(x);
   sz_y = size(y);
   
   %% check to see that the data arrays have the same size
   if (n_x ~= n_y)
     disp(sprintf(['wcsxy2sph.m: input pixel lists do not have same ' ...
                   'length!']));
     return
   end

   %% WE WOULD BRANCH BY MAP PROJECTION TYPE IF WE WANTED MORE THAN ONE HERE
   %% THIS IS FOR MAP TYPE 'TAN'
   theta = pi2 .* ones(size(x));      %%Default is 90 degrees
   r = sqrt(x.^2 + y.^2);
   g = find(r > 0);
   Ng = numel(g);
   if Ng > 0 
     theta(g) = atan2(radeg, r(g));
   end
   phi = atan2(x,-y);
    
   %% Generate input parameters required for next stage.
   wcsparams.phi = phi .* radeg;
   wcsparams.theta = theta .* radeg;
   wcsparams.theta0 = 90.0;
   
   %% Convert from "native" coordinate system to "standard" coordinate
   %% system.
   wcsparams = wcs_rotate(wcsparams);
   longitude = wcsparams.longitude;
   latitude  = wcsparams.latitude;
      
   %% Make longitude have right sign convention
   temp = find(longitude < 0.d0);
   Nneg = numel(temp);
   if (Nneg > 0) 
     longitude(temp) = longitude(temp) + 3.6d2;
   end
   temp = find(longitude >= 3.6d2-1.d-2);
   Nneg = numel(temp);
   if (Nneg > 0) 
     longitude(temp) = longitude(temp) - 3.6d2;
   end

   %% and that's it
   return

end
