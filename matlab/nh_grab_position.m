function nh_grab_position()

  paths = get_paths();

  datafiles = dir(sprintf('%s*.mat',paths.datadir));
  darkfiles = dir(sprintf('%s*.mat',paths.darkdir));
  datafiles = [darkfiles;datafiles];

  nfiles = size(datafiles,1);
  
  darkflag = zeros(nfiles,1);
  darkflag(1:numel(darkfiles)) = 1;
  
  mytime = zeros(nfiles,1);
  mydate = cellstr('');
  myfieldnum = zeros(nfiles,1);
  myfieldflag = zeros(nfiles,1);
  mysunx = zeros(nfiles,1);
  mysuny = zeros(nfiles,1);
  mysunz = zeros(nfiles,1);
  myee = NaN.*zeros(nfiles,1);
  myelong = NaN.*zeros(nfiles,1);
  myelat = NaN.*zeros(nfiles,1);
  
  goodfiles = [3,5,6,7];
  clearedfiles = zeros(4,1);
  
  for ifile=1:nfiles

    clflag = 0;
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));

    if darkflag(ifile) == 1
      
      load(sprintf('%s%s',paths.darkdir,datafiles(ifile).name));
      
    else
      
      load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));

    end
      
    mytime(ifile) = data.header.date_jd-data.header.launch_jd;

    mydate{ifile} = data.header.date_cal;

    myfieldnum(ifile) = -1;
    
    if isfield(data.header,'fieldnum')
      whpl = data.header.fieldnum == goodfiles;
      myfieldnum(ifile) = data.header.fieldnum;
      if clearedfiles(whpl) == 1
	clflag = 1;
      end
      if sum(whpl) > 0 & clflag == 0
	myfieldflag(ifile) = 1;
	clearedfiles(whpl) = 1;
      end
    end
      
    if darkflag(ifile) == 0;
      
      mysunx(ifile) = data.astrom.spcsscx ./ 1.496e8;
      
      mysuny(ifile) = data.astrom.spcsscy ./ 1.496e8;
    
      mysunz(ifile) = data.astrom.spcsscz ./ 1.496e8;

      myee(ifile) = data.coords.sol_elon;
    
      myelong(ifile) = data.coords.ecliptic(1);
      
      myelat(ifile) = data.coords.ecliptic(2);
    
    else

      mysunx(ifile) = data.astrometry.spcsscx ./ 1.496e8;
      
      mysuny(ifile) = data.astrometry.spcsscy ./ 1.496e8;
      
      mysunz(ifile) = data.astrometry.spcsscz ./ 1.496e8;
    
    end
    
  end

  %mytime = [0;mytime];
  %mysunx = [1;mysunx];
  %mysuny = [0;mysuny];
  %mysunz = [0;mysunz];
  
  figure(1); clf
  plot(mytime,mysunx,'o')
  newtime = [0:3500];
  newsunx = interp1(mytime,mysunx,newtime);
  hold on
  plot(newtime,newsunx)
  
  figure(2); clf
  plot(mytime,mysuny,'o')
  newtime = [0:3500];
  newsuny = interp1(mytime,mysuny,newtime);
  hold on
  plot(newtime,newsuny)
  
  figure(3); clf
  plot(mytime,mysunz,'o')
  newtime = [0:3500];
  newsunz = interp1(mytime,mysunz,newtime);
  hold on
  plot(newtime,newsunz)
  
  figure(4); clf
  plot(newsunx,newsuny)
  hold on
  plot(mysunx,mysuny,'o')
  xlim([-30,30])
  ylim([-30,30])
  
  figure(5); clf
  plot(sqrt(newsunx.^2+newsuny.^2),newsunz)
  hold on
  plot(sqrt(mysunx.^2+mysuny.^2),mysunz)
  
  fieldnum = myfieldnum;
  fieldflag = myfieldflag;
  t = mytime;
  x = mysunx;
  y = mysuny;
  r = sqrt(mysuny.^2 + mysunx.^2);
  z = mysunz;
  se = myee;
  elong = myelong;
  elat = myelat;
  
  save('../scratch/nh_position.mat','fieldnum','fieldflag',...
      't','x','y','r','z','se','elong','elat');

  dbstop
  
% end

