% nh_dark_analysis_combined_data()
%shaina thayer July 2019
%anna dignan July 2019

close all;

%Clear temporary variables
clearvars idx;

%Import paths for data location
paths = get_paths_old();
npaths = get_paths_new();

fprintf('Parsing dark files.\n');

%Load all data directories

%dark data sourced from new data directory
ndarkfiles = dir(sprintf('%s*.mat',npaths.darkdir));
%old light data
lightfiles = dir(sprintf('%s*.mat',paths.datadir));
%new light data
nlightfiles = dir(sprintf('%s*.mat',npaths.datadir));

%Preallocate space for variables dark
darktemp = zeros(size(ndarkfiles,1),1);
darkdate = zeros(size(ndarkfiles,1),1);
darksig = zeros(size(ndarkfiles,1),2);
darksig_mean = zeros(size(ndarkfiles,1),2);
darkref = zeros(size(ndarkfiles,1),2);
darkref_med = zeros(size(ndarkfiles,1),2);
darkref_sig = zeros(size(ndarkfiles,1),2);
darkexp = zeros(size(ndarkfiles,1),1);
darkfield = zeros(size(ndarkfiles,1),1);

%Determine which light files correspond to 'good' fields
isgoodold = zeros(numel(lightfiles),1);
isgoodnew = zeros(numel(nlightfiles),1);
%These field numbers were previously determined by going through all the
%files
oldgoodfields = [3,5,6,7];
newgoodfields = [5,6,7,8];

%Check for old light files
for ifile=1:numel(lightfiles)
    load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
    if sum(data.header.fieldnum == oldgoodfields)
        isgoodold(ifile) = 1;
    end
end

%Check for new light files
for ifile=1:numel(nlightfiles)
    load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
    if sum(data.header.fieldnum == newgoodfields)
        isgoodnew(ifile) = 1;
    end
end

%Number of old and new light files corresponding to good fields
numoldlightfiles = sum(isgoodold);
numnewlightfiles = sum(isgoodnew);

%Preallocate space for variables light
lighttemp = zeros((numoldlightfiles+numnewlightfiles),1);
lightfpubtemp = zeros((numoldlightfiles+numnewlightfiles),1);
lightdate = zeros((numoldlightfiles+numnewlightfiles),1);
lightsig = zeros((numoldlightfiles+numnewlightfiles),2);
lightref = zeros((numoldlightfiles+numnewlightfiles),2);
lightref_med = zeros((numoldlightfiles+numnewlightfiles),2);
lightref_sig = zeros((numoldlightfiles+numnewlightfiles),2);
lightexp = zeros((numoldlightfiles+numnewlightfiles),1);
lightfield = zeros((numoldlightfiles+numnewlightfiles),1);
lightlIl = zeros((numoldlightfiles+numnewlightfiles),2);
lightbad = zeros((numoldlightfiles+numnewlightfiles),1);
photcurr = zeros((numoldlightfiles+numnewlightfiles),1);
gall = zeros((numoldlightfiles+numnewlightfiles),1);
galb = zeros((numoldlightfiles+numnewlightfiles),1);
helidist = zeros((numoldlightfiles+numnewlightfiles),1);
%Preallocate space for variables light colorcoded
lightdatenew = zeros((numoldlightfiles+numnewlightfiles),1);
lightdateold = zeros((numoldlightfiles+numnewlightfiles),1);
lighttempnew = zeros((numoldlightfiles+numnewlightfiles),1);
lighttempold = zeros((numoldlightfiles+numnewlightfiles),1);
%Preallocate space for variables light colorcoded dark current
darkcurrnew = zeros((numoldlightfiles+numnewlightfiles),1);
darkcurrold = zeros((numoldlightfiles+numnewlightfiles),1);

%For dark data files
for ifile=1:size(ndarkfiles)
    %Load data and save values
    load(sprintf('%s%s',npaths.darkdir,ndarkfiles(ifile).name));
    
    darktemp(ifile,1) = data.header.ccdtemp;
    darkdate(ifile,1) = data.header.date_jd - data.header.launch_jd;
    darksig(ifile,1) = median(data.dark(:));
    darksig_mean(ifile,1) = mean(data.dark(:));
    darksig(ifile,2) = std(data.dark(:));
    darkref(ifile,1) = mean(data.ref.line);
    darkref_med(ifile,1) = median(data.ref.line);
    darkref_sig(ifile,1) = nh_sigclip(data.ref.line);
    darkref(ifile,2) = std(data.ref.line);
    darkexp(ifile,1) = data.header.exptime;
    %If date in certain range, assign field number
    if darkdate(ifile,1) < 94
        darkfield(ifile,1) = 1;
    end
    if darkdate(ifile,1) > 94 && darkdate(ifile,1) < 96
        darkfield(ifile,1) = 2;
    end
    if darkdate(ifile,1) > 102 && darkdate(ifile) < 103
        darkfield(ifile,1) = 3;
    end
    if darkdate(ifile,1) > 103 && darkdate(ifile,1) < 104
        darkfield(ifile,1) = 4;
    end
    
end

fprintf('Parsing light files.\n');

%For old data files
jfile = 1;
for ifile=1:numel(lightfiles)
    %If file is for a good field, load and save values
    if isgoodold(ifile) == 1
        load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
    
        lighttemp(jfile,1) = data.header.ccdtemp;
        lightfpubtemp(jfile,1) = data.astrom.fpubtemp;
        lightdate(jfile,1) = data.header.date_jd - data.header.launch_jd;
        lightsig(jfile,1) = data.ref.engmean;
        lightsig(jfile,2) = sqrt(data.ref.engmean);
        lightref(jfile,1) = mean(data.ref.line);
        lightref_med(jfile,1) = median(data.ref.line);
        lightref_sig(jfile,1) = nh_sigclip(data.ref.line);
        lightref(jfile,2) = std(data.ref.line);
        lightexp(jfile,1) = data.header.exptime;
        lightfield(jfile,1) = data.header.fieldnum;
        lightlIl(jfile,1) = data.stats.maskmean./data.header.exptime;
        lightlIl(jfile,2) = data.stats.maskstd;
        
        maskim = data.image.calimage.*~data.mask.onemask;
        maskim(maskim==0) = NaN;
        photcurr(jfile,1) = nanmean(nanmean(maskim));
        [l,b]=RA_Dec_to_Gal(data.astrom.spcbrra,data.astrom.spcbrdec);
        gall(jfile,1) = l;
        galb(jfile,1) = b;
        helidist(jfile+numoldlightfiles,1) = sqrt((data.astrom.spcsscx)^2 + (data.astrom.spcsscy)^2 + (data.astrom.spcsscz)^2);

        
        jfile = jfile + 1;
    end
    
darkcurrlight = 2.545.*10.*(1./22).*1e4.*122.*(lighttemp+273).^3.*exp(-6400./(lighttemp+273));

end

%For new data files
jfile = 1;
for ifile=1:numel(nlightfiles)
    %If file is for a good field, load and save values
    if isgoodnew(ifile) == 1
        load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));

        lighttemp(jfile+numoldlightfiles,1) = data.header.ccdtemp;
        lightfpubtemp(jfile+numoldlightfiles,1) = data.astrom.fpubtemp;
        lightdate(jfile+numoldlightfiles,1) = data.header.date_jd - data.header.launch_jd;
        lightsig(jfile+numoldlightfiles,1) = data.ref.engmean;
        lightsig(jfile+numoldlightfiles,2) = sqrt(data.ref.engmean);
        lightref(jfile+numoldlightfiles,1) = mean(data.ref.line);
        lightref_med(jfile+numoldlightfiles,1) = median(data.ref.line);
        lightref_sig(jfile+numoldlightfiles,1) = nh_sigclip(data.ref.line);
        lightref(jfile+numoldlightfiles,2) = std(data.ref.line);
        lightexp(jfile+numoldlightfiles,1) = data.header.exptime;
        lightfield(jfile+numoldlightfiles,1) = data.header.fieldnum;
        lightlIl(jfile+numoldlightfiles,1) = data.stats.maskmean./data.header.exptime;
        lightlIl(jfile+numoldlightfiles,2) = data.stats.maskstd;
        lightbad(jfile+numoldlightfiles,1) = data.header.bad;
        
        maskim = data.image.calimage.*~data.mask.onemask;
        maskim(maskim==0) = NaN;
        photcurr(jfile+numoldlightfiles,1) = nanmean(nanmean(maskim));
        [l,b]=RA_Dec_to_Gal(data.astrom.spcbrra,data.astrom.spcbrdec);
        gall(jfile+numoldlightfiles,1) = l;
        galb(jfile+numoldlightfiles,1) = b;
        helidist(jfile+numoldlightfiles,1) = sqrt((data.astrom.spcsscx)^2 + (data.astrom.spcsscy)^2 + (data.astrom.spcsscz)^2);
        
        jfile = jfile + 1;
    end
darkcurrlight = 2.545.*10.*(1./22).*1e4.*122.*(lighttemp+273).^3.*exp(-6400./(lighttemp+273));

% end
end

%Prepare to save mean values for like fields
nfieldsdark = 4; %Dark data previously divided into 4 fields by date
darktempm = zeros(nfieldsdark,1);
darkerrm = zeros(nfieldsdark,1);
darkdatem = zeros(nfieldsdark,1);
darkrefm = zeros(nfieldsdark,2);

%For each field number, check if dark field matches that number and take
%mean of data values for only that field
for jfield=1:nfieldsdark
    whpl = darkfield == jfield;
    darktempm(jfield) = mean(darktemp(whpl));
    darkerrm(jfield) = std(darktemp(whpl));
    darkdatem(jfield) = mean(darkdate(whpl));
    darkrefm(jfield,1) = sum(darkref(whpl,1) ./ darkref(whpl,2).^2) ./ ...
        sum(1./darkref(whpl,2).^2);
    darkrefm(jfield,2) = sqrt(1./256+std(darkref(whpl,1)).^2);%sqrt(1./sum(1./darkref(whpl,2).^2));
end

%Prepare to save mean values for like fields
nfieldslight = length(oldgoodfields)+length(newgoodfields); %Light data currently comes from 8 different fields
lighttempm = zeros(nfieldslight,1);
lighterrm = zeros(nfieldslight,1);
lightdatem = zeros(nfieldslight,1);
lightrefm = zeros(nfieldslight,2);
lightlIlm = zeros(nfieldslight,2);

%For each old field number, check if light field matches that number and
%take mean of data values for only that field
for jfield=1:numel(oldgoodfields)
    whpl = lightfield(1:numoldlightfiles) == oldgoodfields(jfield);
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

%Create logical vector of zeros to pad for the length of the number of
%good, old light fields
extra = logical(zeros(numoldlightfiles,1));
%For each new field number, take mean of only values corresponding to that
%field 
for jfield=1:numel(newgoodfields)
    whpl = lightfield(numoldlightfiles+1:numnewlightfiles) == newgoodfields(jfield);
    whpl = vertcat(extra,whpl);
    lighttempm(jfield+numel(oldgoodfields)) = mean(lighttemp(whpl));
    lighterrm(jfield+numel(oldgoodfields)) = std(lighttemp(whpl));
    lightdatem(jfield+numel(oldgoodfields)) = mean(lightdate(whpl));
    lightrefm(jfield+numel(oldgoodfields),1) = sum(lightref(whpl,1) ./ lightref(whpl,2).^2) ./ ...
        sum(1./lightref(whpl,2).^2);
    lightrefm(jfield+numel(oldgoodfields),2) = sqrt(1./256 + std(lightref(whpl,1)).^2);%sqrt(1./sum(1./lightref(whpl,2).^2));
    lightlIlm(jfield+numel(oldgoodfields),1) = sum(lightlIl(whpl,1) ./ lightlIl(whpl,2).^2) ./ ...
        sum(1./lightlIl(whpl,2).^2);
    lightlIlm(jfield+numel(oldgoodfields),2) = std(lightlIl(whpl,1));
end

%Ignore any data that are nan
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

%Calculate dark current for temperatures we already have - in e/s/pix
% New gain value from Lauer is 19.4 e/DN
darkcurrlight = 2.2.*2.545.*10.*(1./22).*1e4.*122.*(lighttemp+273.15).^3.*exp(-6400./(lighttemp+273.15));
darkcurrdark = 2.2.*2.545.*10.*(1./22).*1e4.*122.*(darktemp+273.15).^3.*exp(-6400./(darktemp+273.15));
%Load and save dark current
% for ifile=1:numoldlightfiles
%     load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name),'data');
%     data.ref.darkcurr=darkcurrlight(ifile);
%     save(sprintf('%s%s',paths.datadir,lightfiles(ifile).name),'data');
% end
% for ifile=1:numnewlightfiles
%     load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name),'data');
%     data.ref.darkcurr=darkcurrlight(ifile+numoldlightfiles);
%     save(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name),'data');
% end

%Colorcoding dark current for plotting purposes
for ifile=1:numoldlightfiles+numnewlightfiles
    if ifile > numoldlightfiles
       lightdatenew(ifile,1)=lightdate(ifile,1);
       lighttempnew(ifile,1)=lighttemp(ifile,1)+273.15;

    elseif ifile <= numoldlightfiles
       lightdateold(ifile,1)=lightdate(ifile,1);
       lighttempold(ifile,1)=lighttemp(ifile,1)+273.15;

    end
end
for ifile=1:numoldlightfiles+numnewlightfiles
    if ifile > numoldlightfiles
       darkcurrnew(ifile,1)=darkcurrlight(ifile,1);
    elseif ifile <= numoldlightfiles
       darkcurrold(ifile,1)=darkcurrlight(ifile,1);
    end
end

% Convert dark current in e/s/pix to DN/s/pix
darkcurrdn = darkcurrnew/22; % New gain value from Lauer is 19.4 e/DN
% Convert dark current in DN/s/pix to nW/m2/sr
darkcurrnw = darkcurrdn*data.cal.sbconv;
% Calculate mean dark current (few new files) in nW/m2/sr
meandc = mean(darkcurrnw);
% Calculate dark curr err as 20% of dark curr
dc_err = 0.2*meandc;

%Plot photocurrent over mission time for all light data, color coded by new
%vs. old data
% figure;
% k = lightdate < 2500;
% scatter(lightdate(k),photcurr(k),80,'MarkerFaceColor',[235,110,52]/255,'MarkerEdgeColor',[235,110,52]/255);
% hold on;
% scatter(lightdate(~k),photcurr(~k),80,'MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410]);
% xlabel('Days from launch','FontSize',20);
% ylabel('Mean Photocurrent (nW m^-^2 sr^-^1)','FontSize',20);

%Plot photocurrent over mission time for all light data, color coded by
%galactic latitude
% figure;
% scatter(lightdate,photcurr,80,abs(galb));
% colorbar;
% xlabel('Days from launch','FontSize',20);
% ylabel('Mean Photocurrent (nW m^-^2 sr^-^1)','FontSize',20);

figcnt = 1;

% Plot mean of masked raw image (light data) vs. mean of reference pixels
% (also light data) - with x=y line overplotted
figure(figcnt); clf
hold on;
% xlabel('Sigma-Clipped Mean of Reference Pixels')
% ylabel('Mean of Raw Image Pixels')
xlabel('Mean of Raw Image Pixels [DN]')
ylabel('Mean of Reference Pixels [DN]')
% r = scatter(lightref(:,1),lightsig(:,1),'r');
goodlightsig = lightsig(lightbad<1,1);
goodlightref = lightref(lightbad<1,1);
goodlightfield = lightfield(lightbad<1,1);
goodlightdate = lightdate(lightbad<1,1);
goodlightfpubtemp = lightfpubtemp(lightbad<1,1);
% r = scatter(goodlightsig,goodlightref,'r'); % All points red
% r = scatter(goodlightsig,goodlightref,[],goodlightfield); % Points color-coded by field number
% r = scatter(goodlightsig,goodlightref,[],goodlightdate); % Points color-coded by date
r = scatter(goodlightsig,goodlightref,[],goodlightfpubtemp); % Points color-coded by fpub temp
g = colorbar;
% Plot x = y
x = [535:550];
y = [535:550];
xy = plot(x,y,'g');
% xlim([80,4000])
% ylim([535,550]);
% b = scatter(darkref(:,1),darksig_mean(:,1),'b');
% Basic linear fit
[fitobject,gof,output] = fit(goodlightsig,goodlightref,'poly1');
xfit=linspace(535,550);
yfit=(fitobject.p1*xfit + fitobject.p2);
% fitwut = plot(xfit,yfit,'b');
% Fit with slope = 1
mdl = fitlm(goodlightsig,goodlightref-1*goodlightsig,'constant');
yfit=(1*xfit + mean(mdl.Fitted));
% fitwut_constslope = plot(xfit,yfit,'m');
% Calculate rejected points
k = goodlightsig >= goodlightref; %get where goodlightref (y on plot) are under y=x line, (1*goodlightsig-0) is the full math for future ref
goodlightsig_rejection = goodlightsig(k);
goodlightref_rejection = goodlightref(k);
% Fit with slope = 1 and points above x=y rejected
mdl_rejection = fitlm(goodlightsig_rejection,goodlightref_rejection-1*goodlightsig_rejection,'constant');
yfit=(1*xfit + mean(mdl_rejection.Fitted));
fitwut_constslope_rejection = plot(xfit,yfit,'k');
% Robust fit w/ default weighting - this is the good one being used for ref corr
mdl_robust = fitlm(goodlightsig,goodlightref,'RobustOpts','on');
yfit=(mdl_robust.Coefficients{2,1}*xfit + mdl_robust.Coefficients{1,1});
fitwut_robust = plot(xfit,yfit,'c');
% Huber weighting
mdl_robust_huber = fitlm(goodlightsig,goodlightref,'RobustOpts','huber');
yfit=(mdl_robust_huber.Coefficients{2,1}*xfit + mdl_robust_huber.Coefficients{1,1});
fitwut_robust_huber = plot(xfit,yfit,'b');
% legend([xy fitwut fitwut_constslope fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
%     {'X = Y',sprintf('Linear Fit: y = %.3fx + %.3f',fitobject.p1,fitobject.p2), ...
%     sprintf('Linear Fit: y = 1x + %.3f',mean(mdl.Fitted))...
%     ,sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
%     sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
%     sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
%     ,'Location','northeast');
legend([xy fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
    {'X = Y', ...
    sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
    sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
    sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
    ,'Location','northeast');
figcnt = figcnt + 1;

%Plot photocurrent over mission time for all light data, color coded by
%heliocentric distance
figure(figcnt);
scatter(lightdate,photcurr,80,abs(helidist));
colorbar;
xlabel('Days from launch','FontSize',20);
ylabel('Mean Photocurrent (nW m^-^2 sr^-^1)','FontSize',20);
figcnt = figcnt + 1;

%Plot dark current and CCD temp over time, color coded by new vs. old data
fig = figure(figcnt); clf
left_color = [0 0 0];
right_color = [0 0 0];
set(fig, 'defaultAxesColorOrder',[left_color;right_color]);
yyaxis left
p1 = semilogy(lightdatenew(lightdatenew~=0),darkcurrnew(darkcurrnew~=0),'b.', 'MarkerSize',15);
hold on;
dd = semilogy(darkdate, darkcurrdark, 'r.', 'MarkerSize',15);
p2 = semilogy(lightdateold(lightdateold~=0),darkcurrold(darkcurrold~=0),'r.', 'MarkerSize',15);
hold on;
xline(3463,'k:');
xlabel('Days from launch');
ylabel('Dark current (e^-/sec/pixel)');
yyaxis right
r1 = scatter(lightdatenew(lightdatenew~=0),lighttempnew(lighttempnew~=0),'MarkerEdgeColor','none');
r2 = scatter(lightdateold(lightdateold~=0),lighttempold(lighttempold~=0),'MarkerEdgeColor','none');
r3 = scatter(darkdate, darktemp+273.15, 'MarkerEdgeColor','none');
ylabel('CCD Temperature (K)');
legend([p2 p1],{'Pre-Pluto encounter','Pluto encounter and beyond'});
figcnt = figcnt + 1;

% figure(2); clf
% plot(lighttemp,darkcurr,'b.');
% xlabel('CCD Temperature (C)');
% ylabel('Dark current');

figure(figcnt); clf
semilogx(lightdate,lighttemp,'r.');
hold on;
%errorbar(lightdate,lighttemp,0.15.*ones(size(lightdate)),'ro');
semilogx(lightdatem,lighttempm,'ko');
%errorbar(lightdatem,lighttempm,2.*sqrt(lighterrm.^2+(0.15.*ones(size(lightdatem))).^2),'ko');
semilogx(darkdate,darktemp,'b.') ;
hold on;
%errorbar(darkdate,darktemp,0.15.*ones(size(darkdate)),'bo');
semilogx(darkdatem,darktempm,'ko');
%errorbar(darkdatem,darktempm,2.*sqrt(darkerrm.^2 + (0.15.*ones(size(darkdatem))).^2),'ko');
xlim([80,4000]);
xlabel('Days from launch');
ylabel('CCD Temperature (°C)');

x = [darkdate;lightdate];
y = [darktemp;lighttemp];
thismean = median(lighttemp);

z = y - thismean;
f = fit(x,z,'exp1','StartPoint',[-54,0.000109]);
mydates = [50:4000];
myfunc = f.a * exp(f.b .* mydates) + thismean;
semilogx(mydates,myfunc,'b');
hold on;
ylim([-85,-45]);

cover = data.header.cover_jd - data.header.launch_jd;
cover = [cover,cover];
plot(cover,ylim,'k:');
xline(3463,'k:');

lighttemp = lighttemp + 273.15;
darktemp = darktemp + 273.15;
myfunc = myfunc + 273.15;
figcnt = figcnt + 1;

%   save('../scratch/nh_dark_analysis_fig1.mat','lightdate','lighttemp',...
%       'darkdate','darktemp','mydates','myfunc','cover');

figure(figcnt); clf
s1 = semilogx(lightdate,lightref(:,1),'ro');
hold on   ;
%errorbar(lightdate,lightref(:,1),lightref(:,2)./sqrt(256),'ro')
s2 = semilogx(lightdatem,lightrefm(:,1),'kh');
%errorbar(lightdatem,lightrefm(:,1),lightrefm(:,2),'kh')
s3 = semilogx(darkdate,darkref(:,1),'bo') ;
%errorbar(darkdate,darkref(:,1),darkref(:,2)./sqrt(256),'bo')
s4 = semilogx(darkdatem,darkrefm(:,1),'kh');
%errorbar(darkdatem,darkrefm(:,1),darkrefm(:,2),'kh')
xlim([80,4000])
ylim([520,575]);
xlabel('Days from launch');
ylabel('Mean of Reference Pixels');

cover = data.header.cover_jd - data.header.launch_jd;
cover = [cover,cover];
s5 = plot(cover,ylim,'k:');
s6 = xline(3463,'k:');
legend([s1 s2 s3 s4 s5 s6],{'Light Images','Light Images Field Mean','Dark Images','Dark Images Field Mean','Cover Ejected','Pluto Encounter'},'Location','northwest');

meanvref = sum(lightrefm(:,1)./lightrefm(:,2).^2)./sum(1./lightrefm(:,2).^2);
sigvref = std(lightref(:,1))./2;
figcnt = figcnt + 1;

figure(figcnt); clf
hold on;
xlabel('Mean of Reference Pixels')
ylabel('Mean of Raw Image Pixels')
r = scatter(lightref(:,1),lightsig(:,1),'r');
hold on;
b = scatter(darkref(:,1),darksig(:,1),'b');

%vreffit = polyfit(darkref(:,1),darksig(:,1),1);
[a_york, b_york, sigma_ayork, sigma_byork] =...
    york_fit(darkref(:,1)',darksig(:,1)',darkref(:,2)',darksig(:,2)');

vref = [meanvref:0.1:meanvref+11];
vlight = b_york .* vref + a_york;
p1 = plot(vref,vlight,'b');
p2 = plot([meanvref,meanvref],ylim,'r--');
plot([meanvref+sigvref,meanvref+sigvref],ylim,'r:');
plot([meanvref-sigvref,meanvref-sigvref],ylim,'r:');


xlabel('Mean of Reference Pixels');
ylabel('Mean of Image Pixels') ;
legend([r b],{'Light Images','Dark Images',},'Location','southeast');
figcnt = figcnt + 1;

figure(figcnt); clf
plot(darktemp,darkref(:,1),'bo');
hold on

[a_york, b_york, sigma_ayork, sigma_byork] =...
    york_fit(darktemp',darkref(:,1)',0.15.*ones(1,numel(darktemp)),...
    darkref(:,2)');

tccd = [-54.5:0.1:-52.0];
%plot(tccd,b_york.*tccd+a_york,'b');
plot(darktempm,darkrefm(:,1),'kh');
%errorbar(darktempm,darkrefm(:,1),darkrefm(:,2),'kh');
for ifield=1:4
    plot([darktempm(ifield)-0.15,darktempm(ifield)+0.15],...
        [darkrefm(ifield,1),darkrefm(ifield,1)],'k');
end

xlabel('CCD Temperature (K)');
ylabel('Mean of Reference Pixels, Cover On');
figcnt = figcnt + 1;

figure(figcnt); clf
tccd = [-85:1:-50];
plot(tccd,b_york.*tccd+a_york,'b');
hold on
plot(darktempm,darkrefm(:,1),'bh');
%errorbar(darktempm,darkrefm(:,1),darkrefm(:,2),'bh');
for ifield=1:4
    plot([darktempm(ifield)-0.15,darktempm(ifield)+0.15],...
        [darkrefm(ifield,1),darkrefm(ifield,1)],'b');
end
plot(lighttempm,lightrefm(:,1),'rh')
%errorbar(lighttempm,lightrefm(:,1),lightrefm(:,2),'rh')
for ifield=1:numel(lighttempm)
    plot([lighttempm(ifield)-0.15,lighttempm(ifield)+0.15],...
        [lightrefm(ifield,1),lightrefm(ifield,1)],'r');
end
plot(xlim,[meanvref,meanvref],'r--');
plot(xlim,[meanvref+sigvref,meanvref+sigvref],'r:');
plot(xlim,[meanvref-sigvref,meanvref-sigvref],'r:');

xlabel('CCD Temperature (°C)');
ylabel('Mean of Reference Pixels');
figcnt = figcnt + 1;

figure(figcnt); clf
tccd = [-85:1:-45];
darkcurrent = 10.*(1./22).*1e4.*122.*(tccd+273).^3.*exp(-6400./(tccd+273));
plot(tccd,darkcurrent+meanvref,'k');
modelone = darkcurrent+meanvref;
hold on;
darkcurrentm = 2.545*10.*(1./22).*1e4.*122.*(darktempm+273).^3.*exp(-6400./(darktempm+273)) + meanvref;
sum((darkcurrentm - darkrefm(:,1)).^2./darkrefm(:,2).^2);
darkcurrenth = 2.545.*10.*(1./22).*1e4.*122.*(tccd+273).^3.*exp(-6400./(tccd+273));
plot(tccd,darkcurrenth+meanvref,'k:');
modeltwo = darkcurrenth+meanvref;
%errorbar(darktempm,darkrefm(:,1),darkrefm(:,2),'bh');
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

xlabel('CCD Temperature (°C)');
ylabel('Mean of Reference Pixels');
ylim([537,548]);

tccd = tccd + 273.15;
lighttempm = lighttempm + 273.15;
darktempm = darktempm + 273.15;
figcnt = figcnt + 1;

%
%   save('../scratch/nh_dark_analysis_fig6.mat','meanvref',...
%           'sigvref','tccd','modelone','modeltwo','lighttempm',...
% 	  'lightrefm','darktempm','darkrefm');

figure(figcnt); clf
plot(lightrefm(:,1),lightlIlm(:,1),'ro');
%errorbar(lightrefm(:,1),lightlIlm(:,1),lightlIlm(:,2),'ro');
hold on;
for ifield=1:numel(lighttempm)
    plot([lightrefm(ifield,1)-lightrefm(ifield,2),...
        lightrefm(ifield,1)+lightrefm(ifield,2)],...
        [lightlIlm(ifield,1),lightlIlm(ifield,1)],'r');
end

ylabel('Mean of Masked Flight Image (DN)')
xlabel('Mean of Reference Pixels (DN)')

plot(lighttemp,lightref(:,1),'ro');
hold on;

fitx = [darkdate;lightdate];
fity = [darktemp;lighttemp];
thismean = median(lighttemp);
fity = fity - thismean;
f = fit(fitx,fity,'exp1','StartPoint',[219,-0.000036]);
mydates = [80:4000];
myfunc = f.a * exp(f.b .* mydates) + thismean;
semilogx(mydates,myfunc,'b');
ylim([-85,570]);

cover = data.header.cover_jd - data.header.launch_jd;
cover = [cover,cover];
plot(cover,ylim,'k:');

fprintf('done')

