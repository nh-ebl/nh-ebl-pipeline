function nh_dark_analysis()

  mypaths = get_paths_new();

  fprintf('Parsing dark files.\n')
  
  darkfiles = dir(sprintf('%s*.mat','/data/symons/nh_data/dark/'));  
  
  ndarkfiles = numel(darkfiles);

  darktemp = zeros(ndarkfiles,1);
  darkdate = zeros(ndarkfiles,1);
  darksig = zeros(ndarkfiles,2);
  darkref = zeros(ndarkfiles,2);
  darkexp = zeros(ndarkfiles,1);
  darkfield = zeros(ndarkfiles,1);
  
  for ifile=1:ndarkfiles
    
    load(sprintf('%s%s',mypaths.darkdir,darkfiles(ifile).name));

    darktemp(ifile) = data.header.ccdtemp;
    darkdate(ifile) = data.header.date_jd - data.header.launch_jd;
    darksig(ifile,1) = median(data.dark(:));
    darksig(ifile,2) = std(data.dark(:));
    darkref(ifile,1) = mean(data.ref.line);
    darkref(ifile,2) = std(data.ref.line);
    darkexp(ifile) = data.header.exptime;
    if darkdate(ifile) < 94
      darkfield(ifile) = 1;
    end
    if darkdate(ifile) > 94 & darkdate(ifile) < 96
      darkfield(ifile) = 2;
    end
    if darkdate(ifile) > 102 & darkdate(ifile) < 103
      darkfield(ifile) = 3;
    end
    if darkdate(ifile) > 103 & darkdate(ifile) < 104
      darkfield(ifile) = 4;
    end
    
  end

  nfields = 4;
  darktempm = zeros(nfields,1);
  darkerrm = zeros(nfields,1);
  darkdatem = zeros(nfields,1);
  darkrefm = zeros(nfields,2);
    
  for jfield=1:nfields
    whpl = darkfield == jfield;
    darktempm(jfield) = mean(darktemp(whpl));
    darkerrm(jfield) = std(darktemp(whpl));
    darkdatem(jfield) = mean(darkdate(whpl));
    darkrefm(jfield,1) = sum(darkref(whpl,1) ./ darkref(whpl,2).^2) ./ ...
	sum(1./darkref(whpl,2).^2);
    darkrefm(jfield,2) = sqrt(1./256+std(darkref(whpl,1)).^2);%sqrt(1./sum(1./darkref(whpl,2).^2));
  end
    
  disp(sprintf('Parsing light files.'))
  
  lightfiles = dir(sprintf('%s*.mat',mypaths.datadir));  
  
  isgood = zeros(numel(lightfiles),1);
  goodfields = [1,5,6,7,8];
  
  for ifile=1:numel(lightfiles)
    
    load(sprintf('%s%s',mypaths.datadir,lightfiles(ifile).name));
    if sum(data.header.fieldnum == goodfields)
      isgood(ifile) = 1;
    end
  end
  
  lightfiles = sum(isgood);
  %numel(lightfiles);

  lighttemp = zeros(nlightfiles,1);
  lightdate = zeros(nlightfiles,1);
  lightsig = zeros(nlightfiles,2);
  lightref = zeros(nlightfiles,2);
  lightexp = zeros(nlightfiles,1);
  lightfield = zeros(nlightfiles,1);
  lightlIl = zeros(nlightfiles,2);
  
  jfile = 1;
  for ifile=1:numel(lightfiles)
    
    if isgood(ifile) == 1
      load(sprintf('%s%s',mypaths.datadir,lightfiles(ifile).name));

      lighttemp(jfile) = data.header.ccdtemp;
      lightdate(jfile) = data.header.date_jd - data.header.launch_jd;    
      lightsig(jfile,1) = data.ref.engmean;
      lightsig(jfile,2) = sqrt(data.ref.engmean);
      lightref(jfile,1) = mean(data.ref.line);
      lightref(jfile,2) = std(data.ref.line);
      lightexp(jfile) = data.header.exptime;
      lightfield(jfile) = data.header.fieldnum;
      lightlIl(jfile,1) = data.stats.maskmean./data.header.exptime;
      lightlIl(jfile,2) = data.stats.maskstd;
    
      jfile = jfile + 1;
    end
      
  end
  
  nfields = 14;
  lighttempm = zeros(nfields,1);
  lighterrm = zeros(nfields,1);
  lightdatem = zeros(nfields,1);
  lightrefm = zeros(nfields,2);
  lightlIlm = zeros(nfields,2);
    
  for jfield=1:nfields
    whpl = lightfield == jfield;
    lighttempm(jfield) = mean(lighttemp(whpl));
    lighterrm(jfield) = std(lighttemp(whpl));
    lightdatem(jfield) = mean(lightdate(whpl));
    lightrefm(jfield,1) = sum(lightref(whpl,1) ./ lightref(whpl,2).^2) ./ ...
	sum(1./lightref(whpl,2).^2);
    lightrefm(jfield,2) = sqrt(1./256 + std(lightref(whpl,1)).^2);%sqrt(1./sum(1./lightref(whpl,2).^2));
    lightlIlm(jfield,1) = sum(lightlIl(whpl,1) ./ lightlIl(whpl,2).^2) ./ ...
	sum(1./lightlIl(whpl,2).^2);
    lightlIlm(jfield,2) = std(lightlIl(whpl,1));
  end
  
  whpl = ~isnan(lighttempm);
  lighttempm = lighttempm(whpl);
  lighterrm = lighterrm(whpl);
  lightdatem = lightdatem(whpl);
  lightrefmp = lightrefm(whpl,1);
  lightrefmq = lightrefm(whpl,2);
  lightrefm = [lightrefmp,lightrefmq];
  lightlIlmp = lightlIlm(whpl,1);
  lightlIlmq = lightlIlm(whpl,2);
  lightlIlm = [lightlIlmp,lightlIlmq];
  
  fprintf('Making plots.\n')
  
  figure;
  figure(1); clf
  semilogx(lightdate,lighttemp,'r.')
  hold on                           
  %errorbar(lightdate,lighttemp,0.15.*ones(size(lightdate)),'ro')
  %semilogx(lightdatem,lighttempm,'ko');
  %errorbar(lightdatem,lighttempm,2.*sqrt(lighterrm.^2+(0.15.*ones(size(lightdatem))).^2),'ko')
  semilogx(darkdate,darktemp,'b.')  
  %errorbar(darkdate,darktemp,0.15.*ones(size(darkdate)),'bo')
  %semilogx(darkdatem,darktempm,'ko');
  %errorbar(darkdatem,darktempm,2.*sqrt(darkerrm.^2 + (0.15.*ones(size(darkdatem))).^2),'ko');
  xlim([80,4000])
  xlabel('Days from launch')
  ylabel('CCD Temperature (C)')

  x = [darkdate;lightdate];
  y = [darktemp;lighttemp];
  thismean = median(lighttemp);
  z = y - thismean;
  f = fit(x,z,'exp1');
  mydates = [50:4000];  
  myfunc = f.a * exp(f.b .* mydates) + thismean;
  semilogx(mydates,myfunc,'b');
  ylim([-85,-45]);

  cover = data.header.cover_jd - data.header.launch_jd;
  cover = [cover,cover];
  plot(cover,ylim,'k:')
  
  lighttemp = lighttemp + 273.15;
  darktemp = darktemp + 273.15;
  myfunc = myfunc + 273.15;
  
%   save('../scratch/nh_dark_analysis_fig1.mat','lightdate','lighttemp',...
%       'darkdate','darktemp','mydates','myfunc','cover');
  
  figure;
  figure(2); clf
  semilogx(lightdate,lightref(:,1),'ro')
  hold on             
  errorbar(lightdate,lightref(:,1),lightref(:,2)./sqrt(256),'ro')
  semilogx(lightdatem,lightrefm(:,1),'kh');
  errorbar(lightdatem,lightrefm(:,1),lightrefm(:,2),'kh')
  semilogx(darkdate,darkref(:,1),'bo')  
  errorbar(darkdate,darkref(:,1),darkref(:,2)./sqrt(256),'bo')
  semilogx(darkdatem,darkrefm(:,1),'kh');
  errorbar(darkdatem,darkrefm(:,1),darkrefm(:,2),'kh')
  xlim([80,4000])
  xlabel('Days from launch')
  ylabel('Mean of Reference Pixels')
  
  cover = data.header.cover_jd - data.header.launch_jd;
  cover = [cover,cover];
  plot(cover,ylim,'k:')
  
  meanvref = sum(lightrefm(:,1)./lightrefm(:,2).^2)./sum(1./lightrefm(:,2).^2);
  sigvref = std(lightref(:,1))./2;
  
  figure;
  figure(3); clf
  hold on;
  xlabel('ref')
  ylabel('sig')
  scatter(lightref(:,1),lightsig(:,1),'r');
  hold on;
  scatter(darkref(:,1),darksig(:,1),'b');
  
  %vreffit = polyfit(darkref(:,1),darksig(:,1),1);
  [a_york, b_york, sigma_ayork, sigma_byork] =...
      york_fit(darkref(:,1)',darksig(:,1)',darkref(:,2)',darksig(:,2)');
  
  vref = [meanvref:0.1:meanvref+11];
  vlight = b_york .* vref + a_york;
  plot(vref,vlight,'b');
  plot([meanvref,meanvref],ylim,'r--');
  plot([meanvref+sigvref,meanvref+sigvref],ylim,'r:');
  plot([meanvref-sigvref,meanvref-sigvref],ylim,'r:');
  
  
  xlabel('Mean of Reference Pixels')
  ylabel('Mean of Light Pixels')  
  
  figure;
  figure(4); clf
  plot(darktemp,darkref(:,1),'bo')
  hold on
  
  [a_york, b_york, sigma_ayork, sigma_byork] =...
      york_fit(darktemp',darkref(:,1)',0.15.*ones(1,numel(darktemp)),...
      darkref(:,2)');
  
  tccd = [-54.5:0.1:-52.0];
  plot(tccd,b_york.*tccd+a_york,'b')
  plot(darktempm,darkrefm(:,1),'kh')
  errorbar(darktempm,darkrefm(:,1),darkrefm(:,2),'kh')
  for ifield=1:4
    plot([darktempm(ifield)-0.15,darktempm(ifield)+0.15],...
	[darkrefm(ifield,1),darkrefm(ifield,1)],'k');
  end
    
  xlabel('CCD Temperature')
  ylabel('Mean of Reference Pixels, Cover On')
  
  figure;
  figure(5); clf
  tccd = [-85:1:-50];
  plot(tccd,b_york.*tccd+a_york,'b')
  hold on
  plot(darktempm,darkrefm(:,1),'bh')
  errorbar(darktempm,darkrefm(:,1),darkrefm(:,2),'bh')
  for ifield=1:4
    plot([darktempm(ifield)-0.15,darktempm(ifield)+0.15],...
	[darkrefm(ifield,1),darkrefm(ifield,1)],'b');
  end
  plot(lighttempm,lightrefm(:,1),'rh')
  errorbar(lighttempm,lightrefm(:,1),lightrefm(:,2),'rh')
  for ifield=1:numel(lighttempm)
    plot([lighttempm(ifield)-0.15,lighttempm(ifield)+0.15],...
	[lightrefm(ifield,1),lightrefm(ifield,1)],'r');
  end
  plot(xlim,[meanvref,meanvref],'r--');
  plot(xlim,[meanvref+sigvref,meanvref+sigvref],'r:');
  plot(xlim,[meanvref-sigvref,meanvref-sigvref],'r:');

  xlabel('CCD Temperature')
  ylabel('Mean of Reference Pixels')
  
  figure;
  figure(6); clf
  tccd = [-85:1:-45];
  darkcurrent = 10.*(1./22).*1e4.*122.*(tccd+273).^3.*exp(-6400./(tccd+273));
  plot(tccd,darkcurrent+meanvref,'k')
  modelone = darkcurrent+meanvref;
  hold on  
  darkcurrentm = 2.545*10.*(1./22).*1e4.*122.*(darktempm+273).^3.*exp(-6400./(darktempm+273)) + meanvref;
  sum((darkcurrentm - darkrefm(:,1)).^2./darkrefm(:,2).^2)
  darkcurrenth = 2.545.*10.*(1./22).*1e4.*122.*(tccd+273).^3.*exp(-6400./(tccd+273));
  plot(tccd,darkcurrenth+meanvref,'k:')
  modeltwo = darkcurrenth+meanvref;
  errorbar(darktempm,darkrefm(:,1),darkrefm(:,2),'bh')
  for ifield=1:4
    plot([darktempm(ifield)-0.15,darktempm(ifield)+0.15],...
	[darkrefm(ifield,1),darkrefm(ifield,1)],'b');
  end
  
  plot(lighttempm,lightrefm(:,1),'rh')
  errorbar(lighttempm,lightrefm(:,1),lightrefm(:,2),'rh')
  for ifield=1:numel(lighttempm)
    plot([lighttempm(ifield)-0.15,lighttempm(ifield)+0.15],...
	[lightrefm(ifield,1),lightrefm(ifield,1)],'r');
  end
  plot(xlim,[meanvref,meanvref],'r--');
  plot(xlim,[meanvref+sigvref,meanvref+sigvref],'r:');
  plot(xlim,[meanvref-sigvref,meanvref-sigvref],'r:');

  xlabel('CCD Temperature')
  ylabel('Mean of Reference Pixels')
  ylim([537,548])

  tccd = tccd + 273.15;
  lighttempm = lighttempm + 273.15;
  darktempm = darktempm + 273.15;
%   
%   save('../scratch/nh_dark_analysis_fig6.mat','meanvref',...
%           'sigvref','tccd','modelone','modeltwo','lighttempm',...
% 	  'lightrefm','darktempm','darkrefm');

      
      
  figure;    
  figure(7); clf
  plot(lightrefm(:,1),lightlIlm(:,1),'ro')
  hold on
  errorbar(lightrefm(:,1),lightlIlm(:,1),lightlIlm(:,2),'ro')
  for ifield=1:numel(lighttempm)
    plot([lightrefm(ifield,1)-lightrefm(ifield,2),...
	  lightrefm(ifield,1)+lightrefm(ifield,2)],...
	[lightlIlm(ifield,1),lightlIlm(ifield,1)],'r');	
  end
  
  ylabel('Mean of Masked Flight Image (DN)')
  xlabel('Mean of Reference Pixels (DN)')
    
    plot(lighttemp,lightref(:,1),'ro')
  hold on
  
  
  
  
  

  
  fitx = [darkdate;lightdate];
  fity = [darktemp;lighttemp];
  thismean = median(lighttemp);
  fity = fity - thismean;
  f = fit(fitx,fity,'exp1');
  mydates = [80:4000];  
  myfunc = f.a * exp(f.b .* mydates) + thismean;
  semilogx(mydates,myfunc,'b');
  ylim([-85,-45]);

  cover = data.header.cover_jd - data.header.launch_jd;
  cover = [cover,cover];
  plot(cover,ylim,'k:')
  
  
  

%end