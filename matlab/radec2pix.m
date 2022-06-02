function [x,y] = radec2pix(alpha,delta,astrometry,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function radec2pix.m
%%  Sept. 21, 2010
%%  Mike Zemcov
%%
%%  This function takes a CIBER astrometry structure and a set of ra, dec
%%  coordinates and outputs the corresponding pixel numbers. 
%%   
%%  Inputs: astrometry - a CIBER-style astrometry structure
%%          alpha      - ra corresponding to the input pixel vectors
%%          delta      - dec corresponding to the input pixel vectors
%%  Outputs: xpix      - a list of x pixels
%%           ypix      - a list of y pixels
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  p = inputParser;
  
  p.addRequired('alpha',@isnumeric);
  p.addRequired('delta',@isnumeric);
  p.addRequired('astrometry',@isstruct);
  p.addParamValue('cdshift',1,@isnumeric);
  
  p.parse(alpha,delta,astrometry,varargin{:});
  
  alpha      = p.Results.alpha;
  delta      = p.Results.delta;
  astrometry = p.Results.astrometry;
  cdshift    = p.Results.cdshift;
  
%   clear p varargin;
  
   %% define angle constants
%    radeg = 180.0 ./ pi;                 

   %% set up some basic info to pass to the coordinate rotator
   wcsparams.lonpole = astrometry.lonpole;
   wcsparams.latpole = astrometry.latpole;
   wcsparams.crval = [astrometry.crval1,astrometry.crval2];
  
   %% rotate the coordinants
   [xsi,eta] = wcssph2xy(alpha,delta,wcsparams);

   %% make the plate scale matrix
   cd = cdshift.*[[astrometry.cd1_1,astrometry.cd1_2]; ...
                  [astrometry.cd2_1,astrometry.cd2_2]];
   
   %% we need the inverse to go to pixel number
   cdinv = inv(cd);
   
   %% multiply by the relevant vectors
   u = cdinv(1,1) .* xsi + cdinv(1,2) .* eta;
   v = cdinv(2,1) .* xsi + cdinv(2,2) .* eta;

   %% now undo the distortion - note that what is given in the fits
   %% header is only an approximate inversion so this will not be able to
   %% replicate the x,y coords perfectly - for ciber typically the error is
   %% smaller than 1 arcsec over the imager field.  To do this perfectly
   up = u;
   vp = v;
   %% would require an iterative scheme.
   if isfield(astrometry,'ap_order')
        if( astrometry.ap_order == 1 )
            if isfield(astrometry,'ap_0_0')
             up = up + astrometry.('ap_0_0') .* u.^(0) .* v.^(0);
            end
            if isfield(astrometry,'ap_1_0')
             up = up + astrometry.('ap_1_0') .* u.^(1) .* v.^(0);
            end
            if isfield(astrometry,'ap_0_1')
             up = up + astrometry.('ap_0_1') .* u.^(0) .* v.^(1);
            end
            if isfield(astrometry,'ap_1_1')
             up = up + astrometry.('ap_1_1') .* u.^(1) .* v.^(1);
            end
            
            if isfield(astrometry,'bp_0_0')
             vp = vp + astrometry.('bp_0_0') .* u.^(0) .* v.^(0);
            end
            if isfield(astrometry,'bp_1_0')
             vp = vp + astrometry.('bp_1_0') .* u.^(1) .* v.^(0);
            end
            if isfield(astrometry,'bp_0_1')
             vp = vp + astrometry.('bp_0_1') .* u.^(0) .* v.^(1);
            end
            if isfield(astrometry,'bp_1_1')
             vp = vp + astrometry.('bp_1_1') .* u.^(1) .* v.^(1);
            end
        elseif( astrometry.ap_order == 2 )
            if isfield(astrometry,'ap_0_0')
             up = up + astrometry.('ap_0_0') .* u.^(0) .* v.^(0);
            end
            if isfield(astrometry,'ap_1_0')
             up = up + astrometry.('ap_1_0') .* u.^(1) .* v.^(0);
            end
            if isfield(astrometry,'ap_0_1')
             up = up + astrometry.('ap_0_1') .* u.^(0) .* v.^(1);
            end
            if isfield(astrometry,'ap_1_1')
             up = up + astrometry.('ap_1_1') .* u.^(1) .* v.^(1);
            end
            if isfield(astrometry,'ap_0_2')
             up = up + astrometry.('ap_0_2') .* u.^(0) .* v.^(2);
            end
            if isfield(astrometry,'ap_1_2')
             up = up + astrometry.('ap_1_2') .* u.^(1) .* v.^(2);
            end
            if isfield(astrometry,'ap_2_0')
             up = up + astrometry.('ap_2_0') .* u.^(2) .* v.^(0);
            end
            if isfield(astrometry,'ap_2_1')
             up = up + astrometry.('ap_2_1') .* u.^(2) .* v.^(1);
            end
            if isfield(astrometry,'ap_2_2')
             up = up + astrometry.('ap_2_2') .* u.^(2) .* v.^(2);
            end

            if isfield(astrometry,'bp_0_0')
             vp = vp + astrometry.('bp_0_0') .* u.^(0) .* v.^(0);
            end
            if isfield(astrometry,'bp_1_0')
             vp = vp + astrometry.('bp_1_0') .* u.^(1) .* v.^(0);
            end
            if isfield(astrometry,'bp_0_1')
             vp = vp + astrometry.('bp_0_1') .* u.^(0) .* v.^(1);
            end
            if isfield(astrometry,'bp_1_1')
             vp = vp + astrometry.('bp_1_1') .* u.^(1) .* v.^(1);
            end
            if isfield(astrometry,'bp_0_2')
             vp = vp + astrometry.('bp_0_2') .* u.^(0) .* v.^(2);
            end
            if isfield(astrometry,'bp_1_2')
             vp = vp + astrometry.('bp_1_2') .* u.^(1) .* v.^(2);
            end
            if isfield(astrometry,'bp_2_0')
             vp = vp + astrometry.('bp_2_0') .* u.^(2) .* v.^(0);
            end
            if isfield(astrometry,'bp_2_1')
             vp = vp + astrometry.('bp_2_1') .* u.^(2) .* v.^(1);
            end
            if isfield(astrometry,'bp_2_2')
             vp = vp + astrometry.('bp_2_2') .* u.^(2) .* v.^(2);
            end
        else
           [iaX,jaY] = meshgrid(0:astrometry.ap_order,0:astrometry.ap_order);
           thisord_mat = strcat('ap_',arrayfun(@num2str, jaY, 'UniformOutput', 0),'_',arrayfun(@num2str, iaX, 'UniformOutput', 0)); %vectorized strings
           for ia=0:astrometry.ap_order
             for ja=0:astrometry.ap_order
    %            thisord = sprintf('ap_%d_%d',ia,ja);
               thisord = thisord_mat{ia+1,ja+1};
               if isfield(astrometry,thisord)
        %          up = up + getfield(astrometry,thisord) .* u.^(ia) .* v.^(ja);
                 up = up + astrometry.(thisord) .* u.^(ia) .* v.^(ja);
               end
             end
           end
           
    %        vp = v; %already done above
           [iaX,jaY] = meshgrid(0:astrometry.bp_order,0:astrometry.bp_order);
           thisord_mat = strcat('bp_',arrayfun(@num2str, jaY, 'UniformOutput', 0),'_',arrayfun(@num2str, iaX, 'UniformOutput', 0)); %vectorized strings
           for ia=0:astrometry.bp_order
             for ja=0:astrometry.bp_order
    %            thisord = sprintf('bp_%d_%d',ia,ja);
               thisord = thisord_mat{ia+1,ja+1};
               if isfield(astrometry,thisord)
        %          vp = vp + getfield(astrometry,thisord) .* u.^(ia) .* v.^(ja);
                 vp = vp + astrometry.(thisord) .* u.^(ia) .* v.^(ja);
               end
             end
           end
        end
    

   end
    
   %% now add back in the reference pixel - this is backwards because of
   %% how matlab does indexing
   x = vp + astrometry.crpix2; 
   y = up + astrometry.crpix1; 

   return

end
