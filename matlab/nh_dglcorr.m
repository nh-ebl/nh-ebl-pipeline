function nh_dglcorr()

  [iohm,im,iu,id] = textread('lookup/dgl_absorption.txt','%f,%f,%f,%f');
  
  iue = iu - im;
  ide = im - id;
  iohm = iohm - 0.8;
  
  figure(1); clf
  plot(iohm,im,'o');
  hold on
  errorbar(iohm,im,ide,iue,'o');
  
  xp = [0:0.1:20];
  
  fun = @(x,xdata)(x(1).*(1 -exp(x(2).*xdata)));
  
  x0 = [35,-0.15];

  fitcoeffs = zeros(100,2);
  
  for itry=1:100
  
    whpl = iohm > -1;
    yin = im(whpl) + randn(size(im(whpl))) .* (iue(whpl) + ide(whpl)) ./ 2;
    xin = iohm(whpl);
    xerr = randn(size(xin));
    whpu = xerr >= 0;
    xin(whpu) = xin(whpu) + 0.06.*xerr(whpu);
    whpd = xerr < 0;
    xin(whpd) = xin(whpd) + 0.02.*xerr(whpd);

    if 1
      warning('off','all')
      opts = optimset('Display','off');
      x = lsqcurvefit(fun,x0,xin,yin,[],[],opts);
      warning('on','all')
      
      fitcoeffs(itry,1) = x(1);
      fitcoeffs(itry,2) = x(2);
    else
      fo = fitoptions('poly1','Weights',(iue(whpl) + ide(whpl)) ./ 2,...
	  'Lower',[-100 0],'Upper',[100 0]);
      thisfit = fit(xin,yin,'poly1',fo);
      fitcoeffs(itry) = thisfit.p1;
    end
      
  end
  
  whpl = iohm < 10;
  %myline = polyfit(iohm(whpl),im(whpl),1);
  fo = fitoptions('poly1','Weights',(iue(whpl) + ide(whpl)) ./ 2,...
      'Lower',[-100 0],'Upper',[100 0]);
  thisfit = fit(iohm(whpl),im(whpl),'poly1',fo);
  myline = thisfit.p1;
  
  plot(xp,xp.*myline(1));
  
  modmean = mean(fitcoeffs,1);
  modstd = std(fitcoeffs,1);
  
  plot(xp,modmean(1).*(1-exp(modmean(2).*xp)));
  
  % kind of ran out of steam here since none of the fields we ended up using
  % for EBL make it past this cut
  dbstop
  
  mymod = modmean./myline.*ones(size(xp));
  mymoddn = (modmean-modstd)./myline.*ones(size(xp));
  mymodup = (modmean+modstd)./myline.*ones(size(xp));
  
  whpl = xp <10;
  mymod(whpl) = 1;
  mymoddn(whpl) = 1;
  mymodup(whpl) = 1;
  
  figure(2); clf
  plot(xp,mymod)
  hold on
  plot(xp,mymoddn)
  plot(xp,mymodup)
  legend('model','dn','up')
  
  dglcorr.ohm = xp;
  dglcorr.corr = mymod;
  dglcorr.corrdn = mymoddn;
  dglcorr.corrup = mymodup;
  
  save('lookup/nh_dglcorr.mat','dglcorr');


% end

