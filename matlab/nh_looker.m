function nh_looker()

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
  mycrr = zeros(nfiles,1);
  
  goodfiles = [3,5,6,7];%,9,14];
  thissum = 0;
  
  for ifile=1:nfiles

    %disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));

    mydate(ifile) = data.header.date_jd-data.header.launch_jd;
    
    myfieldnum(ifile) = data.header.fieldnum;
    
    mytarget{ifile} = data.header.date_cal;
    
    mysun(ifile) = data.distance.d_sun_NH;
    
    mydt(ifile) = data.distance.d_NH_target;%./1.496e8;;

    myelong(ifile) = data.coords.sol_elon; 
    
    myexp(ifile) = round(data.header.exptime);
    
    myref(ifile) = data.ref.mean;
    
    myeng(ifile) = data.ref.engmean;
    
    myohm(ifile) = data.dgl.ohmmean;
    
    myav(ifile) = 0.4.*data.header.A_V;
    
    mysig(ifile) = data.stats.corrmean;
    
    myisl(ifile) = data.isl.trimean + data.isl.usnowing;

    mydgl(ifile) = data.dgl.dglmean; %- data.dgl.dglerr;
    
    myerr(ifile) = data.stats.correrr;
    
    mycrr(ifile) = data.ref.bias;
    
    mygood = myfieldnum(ifile) == goodfiles;
    
    if sum(mygood) > 0
      isgood(ifile) = 1;
      mystring = sprintf('%f, %d, %s, %f, %f, %f, %d, %f, %f, %f, %f',...
	  mydate(ifile),myfieldnum(ifile),...
	  mytarget{ifile},mysun(ifile),myohm(ifile),myelong(ifile),...
	  mycrr(ifile),...
	  mysig(ifile),myisl(ifile),mydgl(ifile),...
	  mysig(ifile)-myisl(ifile)-mydgl(ifile));
    
      disp(mystring)
      
    end

    %if myfieldnum(ifile) == 7
    %  dbstop
    %end

    
  end

  
  isgood = logical(isgood);
  
  mycorr = zeros(nfiles,1);
  
  figure(1); clf
  plot(mysun(isgood),mysig(isgood)-myisl(isgood)-mydgl(isgood),'o')

  
  mydist = zeros(numel(goodfiles),1);
  mysubmen = zeros(numel(goodfiles),1);
  mysuberr = zeros(numel(goodfiles),1);
  mymean = zeros(numel(goodfiles),1);
  myunc = zeros(numel(goodfiles),1);
  myavp = zeros(numel(goodfiles),1);
  myohmp = zeros(numel(goodfiles),1);
  
  for ifield=1:numel(goodfiles)
    
    whpl = myfieldnum == goodfiles(ifield);
    mydist(ifield) = mean(mysun(whpl));
    thissig = mysig(whpl)-myisl(whpl)-mydgl(whpl);
    thiserr = myerr(whpl);
    mysubmen(ifield) = sum(mysig(whpl)./thiserr.^2)./sum(1./thiserr.^2);
    mysuberr(ifield) = std(mysig(whpl));
    mymean(ifield) = sum(thissig./thiserr.^2)./sum(1./thiserr.^2);
    myunc(ifield) = std(thissig);%./sqrt(sum(whpl)-1);
    myavp(ifield) = sum(myav(whpl)./thiserr.^2)./sum(1./thiserr.^2);
    myohmp(ifield) = sum(myohm(whpl)./thiserr.^2)./sum(1./thiserr.^2);
    
  end
      
  figure(2); clf
  plot(mydist,mymean,'o');
  hold on
  errorbar(mydist,mymean,myunc,'o');
    
  supermean = sum(mymean./myunc.^2)./sum(1./myunc.^2)
  superunc = 1./sqrt(sum(1./myunc.^2))
  
  superav = sum(myavp./myunc.^2)./sum(1./myunc.^2)
  
  supersub = sum(mysubmen./myunc.^2)./sum(1./myunc.^2)
  
  plot(xlim,[supermean,supermean],'b')
  plot(xlim,[supermean+superunc,supermean+superunc],'b:')
  plot(xlim,[supermean-superunc,supermean-superunc],'b:')
  plot(xlim,[0,0],'k:')
  xlabel('Heliocentric Distance (AU)')
  ylabel('Optical EBL (nW/m^2/sr)')
  title('Zemcov et al, Preliminary w/ Stat. Errors')
  

  
  dbstop

% end
