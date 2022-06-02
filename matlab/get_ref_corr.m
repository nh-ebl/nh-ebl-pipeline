% This script independently calculates the appropriate reference correction
% (based on analysis previously done in nh_dark_analysis_combined_data but
% extracted here for clarity and to separate dark current analysis) to be
% used in nh_calibrate

% Can use any or all data sets and compare mean of unmasked raw image data
% to sigma-clipped mean, regular mean, or reported bias level of reference
% pixels and determine reference correction via robust fit

clear all
close all
clc

%Import paths for data location
paths = get_paths_old();
npaths = get_paths_new();
lpaths = get_paths_lauer();
wpaths = get_paths_newest();

%old light data
lightfiles = dir(sprintf('%s*.mat',paths.datadir));
%new light data
nlightfiles = dir(sprintf('%s*.mat',npaths.datadir));
%lauer light data
llightfiles = dir(sprintf('%s*.mat',lpaths.datadir));
%newest light data
wlightfiles = dir(sprintf('%s*.mat',wpaths.datadir));

%Determine which light files correspond to 'good' fields
isgoodold = zeros(numel(lightfiles),1);
isgoodnew = zeros(numel(nlightfiles),1);
isgoodlauer = zeros(numel(llightfiles),1);
isgoodnewest = zeros(numel(wlightfiles),1);

%These field numbers were previously determined by going through all the
%files
oldgoodfields = [3,5,6,7];
% oldgoodfields = []; % Set good fields to none to skip data set
newgoodfields = [5,6,7,8];
% newgoodfields = [];
lauergoodfields = [1,2,3,4,5,6,7];
% lauergoodfields = [];
lauer_exlude_enable = true; %enables skipping of new sequences
newestgoodfields = [2,4,5,6,7,12,15,16,17,19,20,22,23];
% newestgoodfields = [];
newest_exlude_enable = true; %enables skipping of new sequences

%Check for old light files
parpoolobj = gcp('nocreate'); % check for thread pool, which can't use the load call
if isa(parpoolobj,'parallel.ThreadPool')
    delete(parpoolobj); %threads can't use load and will error out
end
parfor ifile=1:numel(lightfiles)
    datatemp = load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == oldgoodfields)
        isgoodold(ifile) = 1;
    end
end

%Check for new light files
parfor ifile=1:numel(nlightfiles)
    datatemp = load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == newgoodfields)
        isgoodnew(ifile) = 1;
    end
end

%Check for lauer light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = 150; %s, time to skip at start of sequence
lauer_exclude = zeros(numel(llightfiles),1); %will fill up with new sequences to ignore
for ifile=1:numel(llightfiles)
    load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
    if sum(data.header.fieldnum == lauergoodfields)
        isgoodlauer(ifile) = 1;
    end
    if lauer_exlude_enable
        if( ifile == 1 )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to skip dynamically
        end
        %skip a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            lauer_exclude(ifile) = 1; %set this to exlude
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            lauer_exclude(ifile) = 1; %set this to exlude
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to skip
        end
    end
end
lauer_exclude = find(lauer_exclude); %get it into indexes

%Check for newest light files
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip_time = 150; %s, time to skip at start of sequence
newest_exclude = zeros(numel(wlightfiles),1); %will fill up with new sequences to ignore
for ifile=1:numel(wlightfiles)
    load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
    if sum(data.header.fieldnum == newestgoodfields)
        isgoodnewest(ifile) = 1;
    end
    if newest_exlude_enable
        if( ifile == 1 )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to skip dynamically
        end
        %skip a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            newest_exclude(ifile) = 1; %set this to exlude
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            newest_exclude(ifile) = 1; %set this to exlude
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to skip
        end
    end
end
newest_exclude = find(newest_exclude); %get it into indexes

%Number of old and new light files corresponding to good fields
numoldlightfiles = sum(isgoodold);
numnewlightfiles = sum(isgoodnew);
numlauerlightfiles = sum(isgoodlauer);
numnewestlightfiles = sum(isgoodnewest);

%Preallocate
lightref = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
lightref_dns = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
lightsig = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
lightsig_dns = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
lightdate = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
lightbad = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);
lightref_mean = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1); 
lightbias = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1); 
lightexptime = zeros((numoldlightfiles+numnewlightfiles+numlauerlightfiles+numnewestlightfiles),1);

%For old data files
fprintf('Loading old data \n');
jfile = 1;
for ifile=1:numel(lightfiles)
    %If file is for a good field, load and save values
    if isgoodold(ifile) == 1
        %Load data files
        load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
        %Save mean of raw data and reference data
        lightref(jfile,1) = nh_sigclip(data.ref.line);
        lightref_dns(jfile,1) = nh_sigclip(data.ref.line)/data.header.exptime;
        lightsig(jfile,1) = data.ref.engmean;
        lightsig_dns(jfile,1) = data.ref.engmean/data.header.exptime;
        lightdate(jfile,1) = data.header.date_jd - data.header.launch_jd;
        lightref_mean(jfile,1) = mean(data.ref.line);
        lightbias(jfile,1) = median(data.ref.line);
        lightexptime(jfile,1) = data.header.exptime;
        jfile = jfile + 1;
    end
end

%For new data files
fprintf('Loading new data \n');
jfile = 1;
for ifile=1:numel(nlightfiles)
    %If file is for a good field, load and save values
    if isgoodnew(ifile) == 1
        %Load data files
        load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
        %Save mean of raw data and reference data
        lightsig(jfile+numoldlightfiles,1) = data.ref.engmean;
        lightsig_dns(jfile+numoldlightfiles,1) = data.ref.engmean/data.header.exptime;
        lightref(jfile+numoldlightfiles,1) = nh_sigclip(data.ref.line);
        lightref_dns(jfile+numoldlightfiles,1) = nh_sigclip(data.ref.line)/data.header.exptime;
        lightdate(jfile+numoldlightfiles,1) = data.header.date_jd - data.header.launch_jd;
        lightbad(jfile+numoldlightfiles,1) = data.header.bad;
        lightref_mean(jfile+numoldlightfiles,1) = mean(data.ref.line);
        lightbias(jfile+numoldlightfiles,1) = data.astrom.biaslevl;
        lightexptime(jfile+numoldlightfiles,1) = data.header.exptime;
        jfile = jfile + 1;
    end
end

%For lauer data files
fprintf('Loading Lauer data \n');
jfile = 1;
for ifile=1:numel(llightfiles)
    %If file is for a good field, load and save values
    if isgoodlauer(ifile) == 1
        %Load data files
        load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
        %Save mean of raw data and reference data
        lightsig(jfile+numoldlightfiles+numnewlightfiles,1) = data.ref.engmean;
        lightsig_dns(jfile+numoldlightfiles+numnewlightfiles,1) = data.ref.engmean/data.header.exptime;
        lightref(jfile+numoldlightfiles+numnewlightfiles,1) = nh_sigclip(data.ref.line);
        lightref_dns(jfile+numoldlightfiles+numnewlightfiles,1) = nh_sigclip(data.ref.line)/data.header.exptime;
        lightdate(jfile+numoldlightfiles+numnewlightfiles,1) = data.header.date_jd - data.header.launch_jd;
        if isfield(data.header,'bad')
            lightbad(jfile+numoldlightfiles+numnewlightfiles,1) = data.header.bad;
        else
            if( any(lauer_exclude == ifile) )
                lightbad(jfile+numoldlightfiles+numnewlightfiles,1) = 1;
            end
        end
        lightref_mean(jfile+numoldlightfiles+numnewlightfiles,1) = mean(data.ref.line);
        lightbias(jfile+numoldlightfiles+numnewlightfiles,1) = data.ref.biaslevl;
        lightexptime(jfile+numoldlightfiles+numnewlightfiles,1) = data.header.exptime;
        jfile = jfile + 1;
    end
end

%For newest data files
fprintf('Loading newest data \n');
jfile = 1;
for ifile=1:numel(wlightfiles)
    %If file is for a good field, load and save values
    if isgoodnewest(ifile) == 1
        %Load data files
        load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
        %Save mean of raw data and reference data
        lightsig(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = data.ref.engmean;
        lightsig_dns(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = data.ref.engmean/data.header.exptime;
        lightref(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = nh_sigclip(data.ref.line);
        lightref_dns(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = nh_sigclip(data.ref.line)/data.header.exptime;
        lightdate(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = data.header.date_jd - data.header.launch_jd;
        if isfield(data.header,'bad')
            lightbad(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = data.header.bad;
        else
            if( any(newest_exclude == ifile) )
                lightbad(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = 1;
            end
        end
        lightref_mean(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = mean(data.ref.line);
        lightbias(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = data.ref.biaslevl;
        lightexptime(jfile+numoldlightfiles+numnewlightfiles+numlauerlightfiles,1) = data.header.exptime;
        jfile = jfile + 1;
    end
end


% Plot mean of masked raw image (light data) vs. mean of reference pixels
% (also light data) - with x=y line overplotted
figure(); clf
hold on;
% Axis labels
xlabel('Mean of Unmasked Raw Image Pixels [DN]')
ylabel('Sigma-Clipped Mean of Reference Pixels [DN]')
% ylabel('Mean of Reference Pixels [DN]')
% ylabel('Median/Robust Mean of Referee Pixels [DN]')
% X data
goodlightsig = lightsig(lightbad<1,1); % Mean of masked raw image
% Y data
goodlightref = lightref(lightbad<1,1); % Sigma-clipped mean of ref pixels
% goodlightref = lightref_mean(lightbad<1,1); % Actual mean
% goodlightref = lightbias(lightbad<1,1); % Recorded bias (median or robust mean)
% Color data
% goodlightfield = lightfield(lightbad<1,1);
% goodlightdate = lightdate(lightbad<1,1);
% goodlightfpubtemp = lightfpubtemp(lightbad<1,1);
% Scatter plot
r = scatter(goodlightsig,goodlightref,20,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30'); % All points red
% r = scatter(goodlightsig,goodlightref,[],goodlightfield); % Points color-coded by field number
% r = scatter(goodlightsig,goodlightref,[],goodlightdate); % Points color-coded by date
% r = scatter(goodlightsig,goodlightref,[],goodlightfpubtemp); % Points color-coded by fpub temp
% g = colorbar; % If needed
% Plot x = y
x = [535:550];
y = [535:550];
xy = plot(x,y,'Color','#7E2F8E');
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
fitwut_constslope_rejection = plot(xfit,yfit,'Color','k');
% Robust fit w/ default weighting - this is the good one being used for ref corr
mdl_robust = fitlm(goodlightsig,goodlightref,'RobustOpts','on');
yfit=(mdl_robust.Coefficients{2,1}*xfit + mdl_robust.Coefficients{1,1});
fitwut_robust = plot(xfit,yfit,'Color','#4DBEEE');
% Huber weighting
mdl_robust_huber = fitlm(goodlightsig,goodlightref,'RobustOpts','huber');
yfit=(mdl_robust_huber.Coefficients{2,1}*xfit + mdl_robust_huber.Coefficients{1,1});
fitwut_robust_huber = plot(xfit,yfit,'Color','#D95319');
% Legend
legend([xy fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
    {'X = Y', ...
    sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
    sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
    sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
    ,'Location','northeast');


%Plot relationship with 538 DN subtracted
figure(); clf
hold on;
% Axis labels
xlabel('Mean of Unmasked Raw Image Pixels - 538 [DN]')
ylabel('Sigma-Clipped Mean of Reference Pixels - 538 [DN]')
% ylabel('Mean of Reference Pixels - 538 [DN]')
% ylabel('Median/Robust Mean of Reference Pixels - 538 [DN]')
% X data
goodlightsig = lightsig(lightbad<1,1)-538; % Mean of Unmasked Raw Image Pixels
% Y data
goodlightref = lightref(lightbad<1,1)-538; % Sigma-clipped mean of ref pix
% goodlightref = lightref_mean(lightbad<1,1)-538; % Actual mean
% goodlightref = lightbias(lightbad<1,1)-538; % Recorded bias (median or robust mean)
% Color data
% goodlightfield = lightfield(lightbad<1,1);
goodlightdate = lightdate(lightbad<1,1);
% goodlightfpubtemp = lightfpubtemp(lightbad<1,1);
% Scatter plot
r = scatter(goodlightsig,goodlightref,20,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30'); % All points red
% r = scatter(goodlightsig,goodlightref,[],goodlightfield); % Points color-coded by field number
% r = scatter(goodlightsig,goodlightref,[],goodlightdate); % Points color-coded by date
% r = scatter(goodlightsig,goodlightref,[],goodlightfpubtemp); % Points color-coded by fpub temp
% g = colorbar; % If needed
% Plot x = y
x = [535:550]-538;
y = [535:550]-538;
xy = plot(x,y,'Color','#7E2F8E');
xlim([0,12])
ylim([-2,12]);
% Basic linear fit
[fitobject,gof,output] = fit(goodlightsig,goodlightref,'poly1');
xfit=linspace(535-538,550-538);
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
fitwut_constslope_rejection = plot(xfit,yfit,'Color','k');
% Robust fit w/ default weighting - this is the good one being used for ref corr
mdl_robust = fitlm(goodlightsig,goodlightref,'RobustOpts','on');
yfit=(mdl_robust.Coefficients{2,1}*xfit + mdl_robust.Coefficients{1,1});
fitwut_robust = plot(xfit,yfit,'Color','#4DBEEE');
% Huber weighting
mdl_robust_huber = fitlm(goodlightsig,goodlightref,'RobustOpts','huber');
yfit=(mdl_robust_huber.Coefficients{2,1}*xfit + mdl_robust_huber.Coefficients{1,1});
fitwut_robust_huber = plot(xfit,yfit,'Color','#D95319');
% Legend
legend([xy  fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
    {'X = Y', ...
    sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
    sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
    sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
    ,'Location','northeast');


%Plot relationship with 538 DN subtracted and normalized for integration
%time - this is the one that gives the ref corr
figure(); clf
hold on;
% Axis labels
xlabel('(Mean of Unmasked Raw Image Pixels - 538)/Exp. Time [DN/s]')
ylabel('(Sigma-Clipped Mean of Reference Pixels - 538)/Exp. Time [DN/s]')
% ylabel('Mean of Reference Pixels - 538 [DN]')
% ylabel('Median/Robust Mean of Reference Pixels - 538 [DN]')
% X data
lightsigsub = lightsig-538; % Mean of Unmasked Raw Image Pixels
% lightsigsub(1:26,1) = lightsigsub(1:26,1)/9.967;
% lightsigsub(27:329,1) = lightsigsub(27:329,1)/9.967;
% lightsigsub(330:649,1) = lightsigsub(330:649,1)/29.9676;
lightsigsub = lightsigsub./lightexptime; %divide by exposure time automagically
goodlightsig = lightsigsub(lightbad<1,1);
% Y data
lightrefsub = lightref-538; % Sigma-clipped mean of ref pix
% lightrefsub(1:26,1) = lightrefsub(1:26,1)/9.967;
% lightrefsub(27:329,1) = lightrefsub(27:329,1)/9.967;
% lightrefsub(330:649,1) = lightrefsub(330:649,1)/29.9676;
lightrefsub = lightrefsub./lightexptime; %divide by exposure time automagically
goodlightref = lightrefsub(lightbad<1,1);
% goodlightref = lightref_mean(lightbad<1,1)-538; % Actual mean
% goodlightref = lightbias(lightbad<1,1)-538; % Recorded bias (median or robust mean)
% Color data
% goodlightfield = lightfield(lightbad<1,1);
goodlightdate = lightdate(lightbad<1,1);
% goodlightfpubtemp = lightfpubtemp(lightbad<1,1);
% Scatter plot
r = scatter(goodlightsig,goodlightref,20,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30'); % All points red
% r = scatter(goodlightsig,goodlightref,[],goodlightfield); % Points color-coded by field number
% r = scatter(goodlightsig,goodlightref,[],goodlightdate); % Points color-coded by date
% r = scatter(goodlightsig,goodlightref,[],goodlightfpubtemp); % Points color-coded by fpub temp
% g = colorbar; % If needed
% Plot x = y
% x = [535:550]-538;
% y = [535:550]-538;
x = [0:0.6];
y = [0:0.6];
xy = plot(x,y,'Color','#7E2F8E');
% xlim([0,0.6])
% ylim([-2,12]);
% Basic linear fit
[fitobject,gof,output] = fit(goodlightsig,goodlightref,'poly1');
xfit=linspace(0,0.6);
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
fitwut_constslope_rejection = plot(xfit,yfit,'Color','k');
% Robust fit w/ default weighting - this is the good one being used for ref corr
mdl_robust = fitlm(goodlightsig,goodlightref,'RobustOpts','on');
yfit=(mdl_robust.Coefficients{2,1}*xfit + mdl_robust.Coefficients{1,1});
fitwut_robust = plot(xfit,yfit,'Color','#4DBEEE');
% Huber weighting
mdl_robust_huber = fitlm(goodlightsig,goodlightref,'RobustOpts','huber');
yfit=(mdl_robust_huber.Coefficients{2,1}*xfit + mdl_robust_huber.Coefficients{1,1});
fitwut_robust_huber = plot(xfit,yfit,'Color','#D95319');
% Legend
legend([xy  fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
    {'X = Y', ...
    sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
    sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
    sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
    ,'Location','northwest');


% Plot mean of masked raw image (light data dns) vs. mean of reference pixels
% (also light data) - with x=y line overplotted - not correct way of
% integrating over exp time
figure(); clf
hold on;
% Axis labels
xlabel('Mean of Unmasked Raw Image Pixels [DN/s]')
ylabel('Sigma-Clipped Mean of Reference Pixels [DN/s]')
% ylabel('Mean of Reference Pixels [DN/s]')
% ylabel('Median/Robust Mean of Reference Pixels [DN/s]')
% X data
goodlightsig = lightsig_dns(lightbad<1,1); % Mean of masked raw image
% Y data
goodlightref = lightref_dns(lightbad<1,1); % Sigma-clipped mean of ref pixels
% goodlightref = lightref_mean(lightbad<1,1); % Actual mean
% goodlightref = lightbias(lightbad<1,1); % Recorded bias (median or robust mean)
% Color data
% goodlightfield = lightfield(lightbad<1,1);
% goodlightdate = lightdate(lightbad<1,1);
% goodlightfpubtemp = lightfpubtemp(lightbad<1,1);
% Scatter plot
r = scatter(goodlightsig,goodlightref,20,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30'); % All points red
% r = scatter(goodlightsig,goodlightref,[],goodlightfield); % Points color-coded by field number
% r = scatter(goodlightsig,goodlightref,[],goodlightdate); % Points color-coded by date
% r = scatter(goodlightsig,goodlightref,[],goodlightfpubtemp); % Points color-coded by fpub temp
% g = colorbar; % If needed
% Plot x = y
x = [0:60];
y = [0:60];
xy = plot(x,y,'Color','#7E2F8E');
% Basic linear fit
[fitobject,gof,output] = fit(goodlightsig,goodlightref,'poly1');
xfit=linspace(0,60);
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
fitwut_constslope_rejection = plot(xfit,yfit,'Color','k');
% Robust fit w/ default weighting - this is the good one being used for ref corr
mdl_robust = fitlm(goodlightsig,goodlightref,'RobustOpts','on');
yfit=(mdl_robust.Coefficients{2,1}*xfit + mdl_robust.Coefficients{1,1});
fitwut_robust = plot(xfit,yfit,'Color','#4DBEEE');
% Huber weighting
mdl_robust_huber = fitlm(goodlightsig,goodlightref,'RobustOpts','huber');
yfit=(mdl_robust_huber.Coefficients{2,1}*xfit + mdl_robust_huber.Coefficients{1,1});
fitwut_robust_huber = plot(xfit,yfit,'Color','#D95319');
% Legend
legend([xy fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
    {'X = Y', ...
    sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
    sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
    sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
    ,'Location','southeast');

%Plot relationship with 538 DN subtracted
figure(); clf
hold on;
% Axis labels
xlabel('Mean of Unmasked Raw Image Pixels - 18 [DN/s]')
ylabel('Sigma-Clipped Mean of Reference Pixels - 18 [DN/s]')
% ylabel('Mean of Reference Pixels - 18 [DN/s]')
% ylabel('Median/Robust Mean of Reference Pixels - 18 [DN/s]')
% X data
goodlightsig = lightsig_dns(lightbad<1,1)-18; % Mean of Unmasked Raw Image Pixels
% Y data
goodlightref = lightref_dns(lightbad<1,1)-18; % Sigma-clipped mean of ref pix
% goodlightref = lightref_mean(lightbad<1,1)-538; % Actual mean
% goodlightref = lightbias(lightbad<1,1)-538; % Recorded bias (median or robust mean)
% Color data
% goodlightfield = lightfield(lightbad<1,1);
goodlightdate = lightdate(lightbad<1,1);
% goodlightfpubtemp = lightfpubtemp(lightbad<1,1);
% Scatter plot
r = scatter(goodlightsig,goodlightref,20,'MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30'); % All points red
% r = scatter(goodlightsig,goodlightref,[],goodlightfield); % Points color-coded by field number
% r = scatter(goodlightsig,goodlightref,[],goodlightdate); % Points color-coded by date
% r = scatter(goodlightsig,goodlightref,[],goodlightfpubtemp); % Points color-coded by fpub temp
% g = colorbar; % If needed
% Plot x = y
x = [0:60]-18;
y = [0:60]-18;
xy = plot(x,y,'Color','#7E2F8E');
xlim([0,60-18])
ylim([0,60-18]);
% Basic linear fit
[fitobject,gof,output] = fit(goodlightsig,goodlightref,'poly1');
xfit=linspace(0-18,60-18);
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
fitwut_constslope_rejection = plot(xfit,yfit,'Color','k');
% Robust fit w/ default weighting - this is the good one being used for ref corr
mdl_robust = fitlm(goodlightsig,goodlightref,'RobustOpts','on');
yfit=(mdl_robust.Coefficients{2,1}*xfit + mdl_robust.Coefficients{1,1});
fitwut_robust = plot(xfit,yfit,'Color','#4DBEEE');
% Huber weighting
mdl_robust_huber = fitlm(goodlightsig,goodlightref,'RobustOpts','huber');
yfit=(mdl_robust_huber.Coefficients{2,1}*xfit + mdl_robust_huber.Coefficients{1,1});
fitwut_robust_huber = plot(xfit,yfit,'Color','#D95319');
% Legend
legend([xy  fitwut_constslope_rejection fitwut_robust fitwut_robust_huber], ...
    {'X = Y', ...
    sprintf('Linear Fit (above y=x rejected): y = 1x + %.3f',mean(mdl_rejection.Fitted)),...
    sprintf('Robust Fit: y = %.3fx + %.3f',mdl_robust.Coefficients{2,1},mdl_robust.Coefficients{1,1}),...
    sprintf('Robust Huber Fit: y = %.3fx + %.3f',mdl_robust_huber.Coefficients{2,1},mdl_robust_huber.Coefficients{1,1})}...
    ,'Location','southeast');

fprint('done');