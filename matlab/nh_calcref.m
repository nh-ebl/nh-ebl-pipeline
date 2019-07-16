function nh_calcref()

  paths = get_paths_new();

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
  mycrr = zeros(nfiles,1);
  
  goodfiles = [1,5,6,7,8];
  
  for ifile=1:nfiles
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    mydate(ifile) = data.header.date_jd;
    
    myfieldnum(ifile) = data.header.fieldnum;
    
    mytarget{ifile} = data.header.target_name;
    
    mysig(ifile) = data.stats.maskmean.*data.header.exptime;
    
    mycrr(ifile) = data.ref.bias;
    
    mygood = myfieldnum(ifile) == goodfiles;
    
    if sum(mygood) == 1
      isgood(ifile) = 1;
        
      mystring = sprintf('%d, %s, %f, %f',myfieldnum(ifile),...
	  mytarget{ifile},mysig(ifile),mycrr(ifile));

      disp(mystring);
    end
    
  end
      
  isgood = logical(isgood);
  
  scatter(mycrr(isgood),mysig(isgood),'r');
  hold on;
  [rho,pval] = corr(mycrr(isgood),mysig(isgood));
  
  [linfit,fiterr] = polyfit(mycrr(isgood),mysig(isgood),1);
  ste = sqrt(diag(inv(fiterr.R)*inv(fiterr.R')).*fiterr.normr.^2./fiterr.df);

  mycorr = linfit(1).*mycrr + linfit(2);

  [rho,pval] = corr(mycrr(isgood),mysig(isgood)-mycorr(isgood));
  
  scatter(mycrr(isgood),mysig(isgood)-mycorr(isgood),'b');
     title('calcref');   
  correction.date = mydate;
  correction.corr = mycorr;
  
  %save('lookup/nh_refcorr.mat','correction');
    
%end