function solf=get_solvefields_params()
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function get_solvefields_params.m
%%  Sept 8, 2010
%%  Mike Zemcov
%%  This function returns stuff we need to know to run solve-fields.
%%  Inputs:
%%  Outputs:
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % these are parameters which don't require any input (ie boolean)
  solf(1).overwrite = '--overwrite'; % overwrite existing files
  solf(1).noplots = '--no-plots'; % supress plot making
  solf(1).sex = '--use-sextractor'; % uses sextractor to find sources,
                                    % slightly more accurate
  solf(1).guesson = '--guess-scale'; % try to guess pixel scale - saves time
  solf(1).verifyoff = '--no-verify'; % don't try to verify wcs before solving - saves time
  
  % these are parameters which require filenames and paths; we set these
  % in the code itself
  solf(1).backend = '--backend-config'; % backend file switch
  solf(1).outdir = '--dir'; % directory to put output files in
  solf(1).outname = '--out'; % output file names
  solf(1).fitsname = '--new-fits'; % name of output fits file

  % these are the angular scalings of the solver
  solf(1).strscalelow = '--scale-low'; % smallest angular scale to
                                        % consider in solver
  solf(1).scalelow = '0.15'; % degrees
  solf(1).strscalehigh = '--scale-high'; % largest angular scale to
                                         % consider in solver
  solf(1).scalehigh = '0.45'; % degrees
  solf(1).strscaleunits = '--scale-units'; % units for scale params
  solf(1).scaleunits = 'degwidth'; % set limits above to degrees
  
  solf(1).depth = '--depth'; % limit depth (number of stars) examined per index
  solf(1).numstars = '20'; % number of stars to consider before moving on

  % these set the nominal pointing and search radius of the solver
  solf(1).ra = '--ra';
  solf(1).dec = '--dec';
  solf(1).radius = '--radius';
  solf(1).defrad = '0.1';  % nominal 0.1 degree search radius - hopefully
                          % we aren't off by more than that!

  % these control the distortion polynomial order 
  solf(1).tweak = '--tweak-order';
  solf(1).tord = '1';
  solf(1).center = '--crpix-center';
  
  % copy same parameters for imager 2
  solf(2) = solf(1);
  
return
