function nh_pluto_looker()

  paths = get_paths();

  datafiles = dir(sprintf('%s*.mat',paths.datadir));

  nfiles = size(datafiles,1);
  
  isgood = zeros(nfiles,1);
  mydate = zeros(nfiles,1);
  myfieldnum = zeros(nfiles,1);
  mytarget = cellstr('');
  mysun = zeros(nfiles,1);
  mydt = zeros(nfiles,1);
  myelong = zeros(nfiles,1);
  myexp = zeros(nfiles,1);
  myref = zeros(nfiles,1);
  myeng = zeros(nfiles,1);
  myohm = zeros(nfiles,1);
  myav = zeros(nfiles,1);
  mysig = zeros(nfiles,1);
  myisl = zeros(nfiles,1);
  mydgl = zeros(nfiles,1);
  myerr = zeros(nfiles,1);
  
  goodfiles = [3,5,6,7,9,14];
  
  for ifile=1:nfiles
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    myfieldnum(ifile) = data.header.fieldnum;
    
    mytarget{ifile} = data.header.target_name;
    
    mysun(ifile) = data.distance.d_sun_NH;
    
    myref(ifile) = data.ref.mean;
    
    myeng(ifile) = data.ref.engmean;
    
    mysig(ifile) = data.stats.calmean;

    myerr(ifile) = data.stats.maskmean;    
    
    
    mystring = sprintf('%d, %s',myfieldnum(ifile),...
	mytarget{ifile});
    

    
    if myfieldnum(ifile)==8
      isgood(ifile) = 1;
    end
    
    disp(mystring)
  end
      
  isgood = logical(isgood);
        
  figure(1); clf
  plot(myerr(isgood),myref(isgood),'bo')
  hold on
  
  [rho,p] = corr(myerr(isgood),myref(isgood))
  
  whpl = strcmp('CHARON',mytarget)
  plot(myerr(whpl),myref(whpl),'ro')
  
  whpl = strcmp('P4',mytarget)
  plot(myerr(whpl),myref(whpl),'ko')
  
  whpl = strcmp('P5',mytarget)
  plot(myerr(whpl),myref(whpl),'mo')
  
  legend('23.3 AU','26.7 AU','29.9 AU')
  
  [thisfit,s] = polyfit(myerr,myref,1)
  
  x = [0.55:0.01:0.85];
  plot(x,x.*thisfit(1)+thisfit(2))
  
  ste = sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df);
  
  dbstop
  
%end