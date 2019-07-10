function nh_darkcruise_analysis()

  % all files are 10s integrations.
  preprocess = 0;
  
  mypaths = get_paths();

  if preprocess

    myfiles = dir(sprintf('%s*_0x633_sci_1.fit',mypaths.darkfitsdir));

    nfiles = numel(myfiles);

    au2km = 1.496e8;
    
    for ifile=1:nfiles
    
      disp(sprintf('On file %d of %d.',ifile,nfiles));
      
      dark_i = ...
	  fitsread(sprintf('%s%s',mypaths.darkfitsdir,myfiles(ifile).name));
    
      info = fitsinfo(sprintf('%s%s',mypaths.darkfitsdir,myfiles(ifile).name));
    
      data.dark = dark_i;      
      data.astrometry = fits_to_astrometry(info);
    
      date=data.astrometry.spcutcjd;
      date_split = strsplit(date,' ');
      data.header.date_jd = str2double(date_split(2));
      data.header.date_cal = data.astrometry.spcutcal;
    
      % extract target names
      data.header.target_name = data.astrometry.spctcb;
  
      data.header.lambdamin = 350e-9;
      data.header.lambdamax = 850e-9;
      data.header.lambda = data.astrometry.pivot ./ 10;
      data.header.launch_jd = 2453755.291667;
      data.header.launch_date = '2006-01-19T19:00:00';
      data.header.cover_jd = 2453976.500000;
      data.header.cover_date = '2006-08-29T00:00:00';
      data.header.ccdtemp = data.astrometry.ccdtemp;  
      data.header.exptime = data.astrometry.exptime;
    
      % extract distance from target to NH
      data.distance.x_NH_target = data.astrometry.spctscx;
      data.distance.y_NH_target = data.astrometry.spctscy;
      data.distance.z_NH_target = data.astrometry.spctscz;
      data.distance.d_NH_target = sqrt(data.distance.x_NH_target.^2 + ...
	  data.distance.y_NH_target.^2 + data.distance.z_NH_target.^2) ./ au2km;
      
      % extract distance from target to sun
      data.distance.x_sun_target = data.astrometry.spctsox;
      data.distance.y_sun_target = data.astrometry.spctsoy;
      data.distance.z_sun_target = data.astrometry.spctsoz;
      data.distance.d_sun_target = sqrt(data.distance.x_sun_target.^2 + ...
	  data.distance.y_sun_target.^2 + ...
	  data.distance.z_sun_target.^2) ./ au2km;
  
      % extract distance from NH to sun
      data.distance.x_sun_NH = data.astrometry.spcsscx;
      data.distance.y_sun_NH = data.astrometry.spcsscy;
      data.distance.z_sun_NH = data.astrometry.spcsscz;
      data.distance.d_sun_NH = sqrt(data.distance.x_sun_NH.^2 + ...
	  data.distance.y_sun_NH.^2 + data.distance.z_sun_NH.^2) ./ au2km;
      
      % strip out the obseravtion time stamp and add it to the header
      timestamp = data.astrometry.spcutcjd(4:end);
      data.header.timestamp = timestamp;
      
      data.header.rawfile = myfiles(ifile).name;
      data.header.met = data.astrometry.met;
      data.header.metstr = sprintf('%010d',data.astrometry.met);
      data.header.inst = data.astrometry.instru;
      data.header.apid = data.astrometry.apid;
      
      data_cal = sprintf('%s%s%s_%s',data.header.date_cal(1:4),...
	  data.header.date_cal(6:7),data.header.date_cal(9:10),...
	  data.header.metstr(1:6));
	  
      datastring = sprintf('%s/%s_%010d_%s_eng_1.fit',data_cal,...
	  data.header.inst,data.header.met,data.header.apid);
      
      data.ref.file = datastring;
      
      ref_i = ...
	  fitsread(sprintf('%sdata/%s',mypaths.engdir,datastring));
      
      data.ref.line = ref_i(:,257);      
      
      % save the resulting data file
%       save(sprintf('%s%s.mat',mypaths.darkdir,timestamp),'data');
    
    end     
  end
     
  myfiles = dir(sprintf('%s*.mat',mypaths.darkdir));  
  
  nfiles = numel(myfiles);

  mytemp = zeros(nfiles,1);
  mydate = zeros(nfiles,1);
  mysig = zeros(nfiles,2);
  myref = zeros(nfiles,2);
  myexp = zeros(nfiles,1);
  
  for ifile=1:nfiles
    
    load(sprintf('%s%s',mypaths.darkdir,myfiles(ifile).name));

    mytemp(ifile) = data.header.ccdtemp;
    mydate(ifile) = data.header.date_jd - data.header.launch_jd;
    mysig(ifile,1) = median(data.dark(:));
    mysig(ifile,2) = std(data.dark(:));
    myref(ifile,1) = mean(data.ref.line);
    myref(ifile,2) = std(data.ref.line);
    myexp(ifile) = data.header.exptime;
    
  end
  
  figure(1); clf
  plot(mydate,myref(:,1));
  
  figure(2); clf 
  plot(myref(:,1),mysig(:,1));
  
  figure(3); clf
  plot(mytemp,myref(:,1));
  
  figure(4); clf
  plot(mydate,myexp);
  
%end
