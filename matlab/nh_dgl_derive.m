function nh_dgl_derive()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %%  part 1: fit for scaling constant
   
  % read in bandpass
  [lorri_lambda,lorri_trans] = textread('lookup/nh_lorri_bandpass.txt','%f %f');
 
  % read in dgl c(lambda) factor from ZDA04 model
  [dgl_lambda,dgl_conv] = textread('lookup/nh_dgl.txt','%f, %f');
  
  dgl_lambda = 1e3 .* dgl_lambda;
  
  lambda = [300:1000];
  
  % interpolate the dgl model to our grid spacing
  myconv = interp1(dgl_lambda,dgl_conv,lambda,'spline');

  % make vectors of measurements from Ienaka 2013
  othermeas = [1.6,2.3,4.0,3.4,2.2,1.4,2.6,2.4,4.3,3.0,...
	2.1,3.3,4.2,2.7,2.1,4.6,1.2,2.3];
  othererrs = 5.81.*[0.1,0.1,0.3,0.2,0.3,1.1,1.6,0.1,0.1,0.1,...
	0.1,0.4,0.7,0.7,0.1,0.1,0.1,0.1];
  % factor 5.81 is required to make chi^2_min = 1.0
  otherlambda = [0.44,0.49,0.55,0.65,0.47,0.45,0.65,0.45,0.65,0.90,...
	0.46,0.53,0.63,0.83,0.44,0.64,0.45,0.65];
  
  % these are just the Ienaka points, rather than all of them
  %othermeas = [1.61,2.25,4.00,3.37];
  %othererrs = 2.7*[0.11,0.14,0.28,0.21];
  %otherlambda = [0.44,0.49,0.55,0.65];
  
  % this scales from b(lambda) to c(lambda)
  othermeas = othermeas .* 100 ./ otherlambda * 1e-3;
  othererrs = othererrs .* 100 ./ otherlambda * 1e-3;
    
  % show the measurements 
  figure(1); clf
  plot(otherlambda,othermeas,'bo')
  hold on
  errorbar(otherlambda,othermeas,othererrs,'bo')
  
  % make a vector of scaling factors to minimize against (best number will
  % be ~0.5)
  scal = [0:0.01:2];
  thischisq = zeros(numel(scal),1);
  % interpolate the dgl model onto the measurement lambda grid for comparison
  chidgl = interp1(dgl_lambda./1000,dgl_conv,otherlambda,'spline');

  % plot the model on the data; without scaling it will be bad
  plot(otherlambda,chidgl,'ro:')

  % number of degrees of freedom in the fit
  ndof = numel(otherlambda)-1;
  % now compute the chi^2 of the model against the data
  for is=1:numel(scal)
    thischisq(is) = sum((scal(is).*chidgl - othermeas).^2./othererrs.^2)./ndof;
  end

  % now find just the points where delta (chi^2) = 1;
  whplscal = scal(thischisq < 2);
  
  % make a vector which picks out the minimum chi^2 value, and then the two
  % delta(chi^2) values.  This is finding the values of c consistent with
  % the data
  scalmd = [scal(thischisq == min(thischisq)),whplscal(1),whplscal(end)];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %%  part 2: fit for scale the ZDA04 model by the appropriate amount
  
  % this was the old way of doing this
  if 0
    whpl = abs(lambda-640) == min(abs(lambda-640));
  
    % factor 0.357 comes from Brandt & Draine
    myconv = myconv ./ myconv(whpl) .* 0.357;
  end
    
  % ok, compute the band-weighted integral of c over the LORRI bandpass;
  % this is cbar
  for ix=1:3
    cbar(ix) = sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0) .* ...
	myconv.*scalmd(ix)) ./ ...
	sum(interp1(lorri_lambda,lorri_trans,lambda,'linear',0.0));
  end
    
  % this stuff is required if you actually wanted an estimate of DGL 
  % versus I_100
  if 0
    myfreq = 3e8 ./ 100e-6;
    thisiris = [0:10]; % MJy/sr
    thisdgl = (thisiris) .* 1e-20 .* myfreq .* myconv .* 1e9; % nW m^-2 sr^-1
  end
  
  % tell me what the values of cbar are so I can write them down for later
  disp(sprintf('cbar = %5.3f,-%5.3f,+%5.3f',cbar(1),cbar(1)-cbar(2),...
      cbar(3)-cbar(1)))
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %%  part 3a: find the supported values of g from Draine 2003
  
  % read in the Sano g vs b measurments
  [draine_invlambda,draine_g] = ...
      textread('lookup/dgl_g.txt','%f, %f');
  
  draine_lambda = 1000./draine_invlambda;
  
  % interpolate to our wavelengths of interest
  g_trans = interp1(lorri_lambda,lorri_trans,draine_lambda,'linear',0.0);
  
  % the allowable g's are those that are band-weighted averages.
  g(1) = sum(g_trans.*draine_g)./sum(g_trans);
  g(2) = g(1) - 0.1;
  g(3) = g(1) + 0.1;

  % print those out for me to see
  disp(sprintf('g = %5.3f,-%5.3f,+%5.3f',g(1),g(1)-g(2),g(3)-g(1)))
  
  % a range of plausible b values
  thisb = [20:90];
  % normalized to 25 degrees to match Lillie + Witt 1976
  whpl = thisb == 25;
  % this is the attenuation factor as a function of g and b.
  thiscorrfunc = (1 - 1.1 * g(1) * sqrt(sin(thisb .* pi ./ 180)));
  
  % tell me what the scaling amplitude is to normalize this function to 
  % 41 degrees.
  disp(sprintf('A = %6.3f',thiscorrfunc(whpl)));

  % normalize to that to show me
  thiscorrfunc = thiscorrfunc ./ thiscorrfunc(whpl);

  % show me
  figure(1); clf
  plot(thisb,thiscorrfunc)
  hold on
  plot(xlim,[1,1],'b:')
  plot([41,41],ylim,'b:')
  
  % and this is just a little test of what that function is for each of the
  %four LORRI fields
  myb = [85.7384,28.4141,57.6861,62.0328];
  mycorr = interp1(thisb,thiscorrfunc,myb);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %%  part 3b: find the supported values of g from Sano
  
  % read in the Sano g vs b measurments
  [sano_lambda,sano_g,sano_gdn,sano_gup] = ...
      textread('lookup/sano_g.txt','%f, %f, %f, %f');
  
  % interpolate to our wavelengths of interest
  g_trans = interp1(lorri_lambda,lorri_trans,sano_lambda,'linear',0.0);
  
  % the allowable g's are those that are band-weighted averages.
  g(1) = sum(g_trans.*sano_g)./sum(g_trans);
  g(2) = sum(g_trans.*sano_gdn)./sum(g_trans);
  g(3) = sum(g_trans.*sano_gup)./sum(g_trans);

  % print those out for me to see
  disp(sprintf('g = %5.3f,-%5.3f,+%5.3f',g(1),g(1)-g(2),g(3)-g(1)))
  
  % a range of plausible b values
  thisb = [20:90];
  % normalized to 41 degrees to match Ienaka 2013
  whpl = thisb == 41;
  % this is the attenuation factor as a function of g and b.
  thiscorrfunc = (1 - 1.1 * g(1) * sqrt(sin(thisb .* pi ./ 180)));
  
  % tell me what the scaling amplitude is to normalize this function to 
  % 41 degrees.
  disp(sprintf('A = %6.3f',thiscorrfunc(whpl)));

  % normalize to that to show me
  thiscorrfunc = thiscorrfunc ./ thiscorrfunc(whpl);

  % show me
  figure(1); clf
  plot(thisb,thiscorrfunc)
  hold on
  plot(xlim,[1,1],'b:')
  plot([41,41],ylim,'b:')
  
  % and this is just a little test of what that function is for each of the
  %four LORRI fields
  myb = [85.7384,28.4141,57.6861,62.0328];
  mycorr = interp1(thisb,thiscorrfunc,myb);
  
  dbstop
   
% end
