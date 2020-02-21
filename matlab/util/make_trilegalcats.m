function make_trilegalcats()

  paths = get_paths_new();

  smallfields = [1,5,6,7,8];
  nfields = 5;
    
  %% need to account for B and R bands to improve 
  %% magnitude accuracy.
  lambda_mag = [438,545,641,798];
  [filt_lambda,filt_trans] = textread('lookup/nh_lorri_bandpass.txt',...
      '%f %f');
  % convert from column to row vector
  filt_lambda = filt_lambda';
  % convert from column to row vector
  filt_trans = filt_trans'./max(filt_trans);
  
  for ifield=1:1%nfields

    disp(sprintf('On field %d of %d.',ifield,nfields));
    
    indir = dir(sprintf('%s%02d/*.dat',paths.tridir,ifield));

    nfiles = numel(indir);
    
    for ifile=1:nfiles
      
      disp(sprintf('On file %d of %d.',ifile,nfiles));
      
      disp(sprintf('Reading file...'));
      
      % step 1: read in the trilegal catalog information 
      filename = sprintf('%s%02d/%s',paths.tridir,ifield,...
	  indir(ifile).name);

      [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,U,B,V,R,I,J,H,K,a11] = ...
	  textread(filename,['%d %f %f %f %f %f %f %f %f %f %f '...
	    '%f %f %f %f %f %f %f %f %f'],'commentstyle','shell');

      clear a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 U J H K a11;
      
      nstars = min([numel(B),numel(V),numel(R),numel(I)]);
      
      B = B(1:nstars);
      V = V(1:nstars);
      R = R(1:nstars);
      I = I(1:nstars);
      
      whpl = ifield == smallfields;
      if sum(whpl) == 0
	nrands = nstars;
	nstars = round(0.0847 .* nstars);
	ind = randperm(nrands,nstars);
	B = B(ind);
	V = V(ind);
	R = R(ind);
	I = I(ind);
      end
	      
      disp(sprintf('Calculating magnitudes...'));
      
      mag_lorri_filt = zeros(nstars,1);
      mag_lorri_center = zeros(nstars,1);
      
      for jmag=1:nstars
	%mags = [B(jmag),V(jmag),R(jmag),I(jmag)];  
	fine_mags_coeffs = polyfit(lambda_mag,...
	    [B(jmag),V(jmag),R(jmag),I(jmag)],2);
	filt_mags = fine_mags_coeffs(1) .* filt_lambda.^2 + ...
	    fine_mags_coeffs(2) .* filt_lambda + fine_mags_coeffs(3);
	mag_lorri_filt(jmag) = sum(filt_mags .* filt_trans) ./ sum(filt_trans);
	mag_lorri_center(jmag) = fine_mags_coeffs(1) .* nh_lambda().^2 + ...
	    fine_mags_coeffs(2) .* nh_lambda() + fine_mags_coeffs(3);

	%figure(1); clf
	%plot(lambda_mag,mags,'o')
	%hold on
	%plot([607.6,nh_lambda()],[mag_lorri_filt(jmag),mag_lorri_center(jmag)],'r+')
	%pause()
      end
      
      Vp(ifile).B = B; 
      Vp(ifile).V = V;
      Vp(ifile).R = R;
      Vp(ifile).I = I;
      Vp(ifile).mlf = mag_lorri_filt;
      Vp(ifile).mlc = mag_lorri_center;
      Vp(ifile).mll = [607.6,nh_lambda()];
      
    end

    V = Vp;
    
    save(sprintf('lookup/isltrilegal_%02d.mat',ifield),'V','-v7.3');
    
  end

% end
