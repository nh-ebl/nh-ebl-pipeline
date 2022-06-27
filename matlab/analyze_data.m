clear all
close all

%get paths for new data files or old data files
% paths = get_paths_old_ghosts();
% paths = get_paths_old();
paths = get_paths_new();
datadir = paths.datadir;

%load directory of saved mat files for all images
datafiles = dir(sprintf('%s*.mat',datadir));

% Preallocate
mydate = zeros(size(datafiles));
mytemp = zeros(size(datafiles));
mymean = zeros(size(datafiles));
myref = zeros(size(datafiles));
myeng = zeros(size(datafiles));
myisl = zeros(size(datafiles));

photcurr = zeros(size(datafiles));
checkmax = zeros(size(datafiles));
checkmin = zeros(size(datafiles));

ghostdist = zeros(size(datafiles));
diffghostcheck = zeros(size(datafiles));
diffghostreal = zeros(size(datafiles));
I1 = zeros(size(datafiles));
I2 = zeros(size(datafiles));
im_hist = zeros(size(datafiles));
im_ref = zeros(size(datafiles));
im_mean = zeros(size(datafiles));
im_mostprob = zeros(size(datafiles));
maskmeanold = zeros(size(datafiles));
maskmean = zeros(size(datafiles));
jailbar_type = zeros(size(datafiles));
exposure_time = zeros(size(datafiles));
field = zeros(size(datafiles));
date = zeros(size(datafiles));
im_num = zeros(size(datafiles));
ccdtemp = zeros(size(datafiles));
fpubtemp = zeros(size(datafiles));
date_jd = zeros(size(datafiles));
bad = zeros(size(datafiles));

%%
%Plot all ghost masks together
% h = figure(1);
% clf;
% im = zeros(256);
% imagesc(im);
% hold on;
% pbaspect([1 1 1]);
% xlabel('LORRI X Pixels');
% ylabel('LORRI Y Pixels');

% Diff ghost area mask
[xgrid, ygrid] = meshgrid(1:256, 1:256);
mask = ((xgrid-(128)).^2 + (ygrid-(128)).^2) <= 58.2.^2;
ghostdiffmask = zeros(256);
ghostdiffmask(mask) = 1;
% masksize = sum(ghostdiffmask,'all');
% outersize = 256*256-masksize;

% subplot_setup = [3,4]; % [# Rows, # Columns] for subplot
% subplot_tot = prod(subplot_setup); % Number of plots per subplot
% subplot_cntr = 0; % Counts for the subplot
% valu_lims = [-250,250]; % Limits the plot dynamic range
% origCmap = parula(4096); % Based on https://www.mathworks.com/matlabcentral/answers/307318-does-matlab-have-a-nonlinear-colormap-how-do-i-make-one#comment_695509
% dataMax = max(valu_lims);
% dataMin = min(valu_lims);
% centerPoint = 0;
% scalingIntensity = 6;
% x = 1:length(origCmap); 
% x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
% x = scalingIntensity * x/max(abs(x));
% x = sign(x).* exp(abs(x));
% x = x - min(x); x = x*2047/max(x)+1; 
% newCmap = interp1(x, origCmap, 1:2048);
% newCmap(end,:) = [1,1,1]; % Make last value in cmap white (yellow end)
% newCmap(1,:) = [1,1,1]; % Make 1st valuein cmap white (blue end)

%%
%loop through all images
for ifile=1:size(datafiles,1)
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    %load individual mat file for image
    load(sprintf('%s%s',datadir,datafiles(ifile).name));
    disp(data.header.rawfile)
    
    %% Diffuse ghost subtraction check
    % Plot new ghost mask - im>0 for all masks = 1, im for weighted by
    % repeated location
%     im = im + data.mask.ghostmask;
%     imagesc(im>0);
    
    % Calculate distance from ghost to center of fov + ghost radius
    ghostdist(ifile) = data.ghost.ghostdistcent + 21.5;
    
    % Calculate check of diffuse ghost sum from data (ghost region - not
    % ghost region)
    
    % Ghost region
%     masksize = sum(data.mask.onemask.*ghostdiffmask,'all'); % Mask-based number of pixels
%     I2(ifile) = 1/masksize*sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all');
    %Preallocate for edge values of peak histogram value
%     minmaxedge = zeros(2,1);
    %Save only ghost area of image
    ghost = data.image.calimage(~data.mask.onemask & ghostdiffmask);
    %Save the bin edges and bin counts for a histogram of the ghost
%     [N,edges] = histcounts(ghost);
    %Find the bin with maximum counts
%     [M,I] = max(N);
    %Calculate edge values for that bin
%     minmaxedge(1,1) = edges(I);
%     minmaxedge(2,1) = edges(I+1);
    %Pixel value is average of those bin edges
    I2(ifile) = mean(ghost);
    
    % Non-ghost region
%     outersize = sum(~data.mask.onemask.*~ghostdiffmask,'all'); % Mask-based number of pixels
%     I1(ifile) = 1/outersize*sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all');
    %Preallocate for edge values of peak histogram value
%     minmaxedge = zeros(2,1);
    %Save only ghost area of image
    ghost = data.image.calimage(~data.mask.onemask & ~ghostdiffmask);
    %Save the bin edges and bin counts for a histogram of the ghost
%     [N,edges] = histcounts(ghost);
    %Find the bin with maximum counts
%     [M,I] = max(N);
    %Calculate edge values for that bin
%     minmaxedge(1,1) = edges(I);
%     minmaxedge(2,1) = edges(I+1);
    %Pixel value is average of those bin edges
    I1(ifile) = mean(ghost);
    
    % Difference between thost and non-ghost region
    diffghostcheck(ifile) = I2(ifile) - I1(ifile);
%     diffghostcheck(ifile) = 1/(256*256)*(1/masksize*sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all')-1/outersize*sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all'));
%     diffghostcheck(ifile) = 1/(256*256)*(sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all')-sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all'));
    
    % Diffuse ghost contribution from model
    diffghostreal(ifile) = data.ghost.diffusesub;
%     diffghostreal(ifile) = data.ghost.diffusesub/masksize;

    % Calculate peak hist of masked image
    %Preallocate for edge values of peak histogram value
%     minmaxedge = zeros(2,1);
    %Save only ghost area of image
%     im = data.image.calimage(~data.mask.onemask);
    uncorr_cal = (((data.data./ data.header.exptime)).* data.cal.sbconv);
    im = uncorr_cal(~data.mask.onemask);
    %Find indices where pixel values are > 0
    %idx = im~=0;
%     nbins_rice = round(2*(length(im(idx))^(1/3))); %rice
%     nbins_sturges = ceil(log2(length(im(idx))))+1; %sturges
%     nbins_sqrt = round(sqrt(length(im(idx)))); %sqrt
%     nbins_scott = round((max(im(idx))-min(im(idx)))/(3.5*(std(im(idx))/length(im(idx))^(1/3)))); %scott
%     nbins_fd = round((max(im(idx))-min(im(idx)))/(2*(iqr(im(idx))/length(im(idx))^(1/3)))); %fd
    %Save the bin edges and bin counts for a histogram of the ghost 
    
    %Plot histogram of masked calibrated image
%     h = figure(1);
%     clf;
%     set(h,'visible','off');
%     g = histogram(im);
%     g.BinWidth = 15;
%     hold on;
%     xl1 = xline(median(im),'r');
%     xl1.LineWidth = 2;
%     xl2 = xline(mean(im),'b');
%     xl2.LineWidth = 2;
%     legend([xl1,xl2],{'Median','Mean'},'Location','northeast');
%     title(sprintf('%s',data.header.rawfile));
%     xlabel('Masked Image Intensity [nW m^{-2} sr^{-1}]');
%     ylabel('N');
%     ext = '.png';
%     imagename = sprintf('%s%s%s',paths.maskhistdir,data.header.timestamp,ext);
%     print(h,imagename, '-dpng');
    
%     [~,I] = max(N); % Find bin with max counts
%     minmaxedge(1,1) = edges(I); % Get the edge values for that bin
%     minmaxedge(2,1) = edges(I+1);
%     idx2 = (im >= edges(I)) & im < (edges(I+1)); % Get the values in the bin edges (uses histcount bin edge logic)
%     nbins_rice = round(2*(length(im(idx&idx2))^(1/3))); %rice
%     nbins_sturges = ceil(log2(length(im(idx&idx2))))+1; %sturges
%     nbins_sqrt = round(sqrt(length(im(idx&idx2)))); %sqrt
%     nbins_scott = round((max(im(idx&idx2))-min(im(idx&idx2)))/(3.5*(std(im(idx&idx2))/length(im(idx&idx2))^(1/3)))); %scott
%     nbins_fd = round((max(im(idx&idx2))-min(im(idx&idx2)))/(2*(iqr(im(idx&idx2))/length(im(idx&idx2))^(1/3)))); %fd
%     
%     [N,edges] = histcounts(im(idx&idx2),nbins_sturges); % Redo histcounts with only edges, let auto alg work
%     
%     edges = linspace(edges(1),edges(end),length(edges)*10); % Boost edge fidelity
% %     bin_width = mean(diff(edges))
%     nbins_rice = round(2*(length(im(idx&idx2))^(1/3))); %rice
%     nbins_sturges = ceil(log2(length(im(idx&idx2))))+1; %sturges
%     nbins_sqrt = round(sqrt(length(im(idx&idx2)))); %sqrt
%     nbins_scott = round((max(im(idx&idx2))-min(im(idx&idx2)))/(3.5*(std(im(idx&idx2))/length(im(idx&idx2))^(1/3)))); %scott
%     nbins_fd = round((max(im(idx&idx2))-min(im(idx&idx2)))/(2*(iqr(im(idx&idx2))/length(im(idx&idx2))^(1/3)))); %fd
    
%     [N,edges] = histcounts(im(idx&idx2),edges); % Redo histcounts again with higher fidelity bin edges based on auto alg
    %Find the bin with maximum counts
%     [M,I] = max(N);
%     %Calculate edge values for that bin
%     minmaxedge(1,1) = edges(I);
%     minmaxedge(2,1) = edges(I+1);
%     %Pixel value is average of those bin edges
    im_hist(ifile) = mean(im);
%     if im_hist(ifile) > 11
%         fprintf('ahhhh');
%     end
    im_ref(ifile) = data.ref.bias;
    im_mean(ifile) = data.stats.corrmean;
    im_mostprob(ifile) = data.stats.corrmostprob;
    maskmeanold(ifile) = data.stats.maskmean_dns_frommask;
    maskmean(ifile) = data.stats.maskmean_dns_fromjb;
    jailbar_type(ifile) = data.header.jailbar_type;
    exposure_time(ifile) = data.header.exptime;

    % Plot masked, calibrated image with red circle indicating ghost region
    % Title includes sum and npix of inner and outer regions, and difference of
    % inner - outer, which estimates diff ghosts
    
    %----- Multi-Image Plot Saver -----
%     if( mod(ifile-1,subplot_tot) == 0 )        
%         h = figure(1);
%         if( ifile == 1 )
%             set(h, 'units', 'normalized'); % Force full screen
%             set(h, 'outerposition', [0 0 1 1]); % Force full screen
% %             set(h,'visible','off'); % Causes figure sizing problems
%         end
%         clf;
%         % tiledlayout instead of subplot per https://www.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html
%         tiledlayout(h,subplot_setup(1),subplot_setup(2),'TileSpacing','Compact','Padding','Compact'); 
%     end
% %     subplot(subplot_setup(1),subplot_setup(2),mod(ifile-1,subplot_tot)+1); % Tiled layout instead
%     nexttile; % Tiled layout instead
%     imagesc(data.image.calimage.*~data.mask.onemask);
%     a = colorbar;
%     a.Label.String = 'Intensity [nW m^{-2} sr^{-1}]';
%     caxis([min(valu_lims),max(valu_lims)]);
% %     colormap(newCmap);
%     hold on;
%     pbaspect([1 1 1]);
%     xlabel('LORRI X Pixels');
%     ylabel('LORRI Y Pixels');
%     th = 0:pi/50:2*pi;
%     xunit = 58.2*cos(th)+128;
%     yunit = 58.2*sin(th)+128;
%     plot(xunit,yunit,'r');
%     title([num2str(ifile),' ',num2str(data.header.timestamp)])
%     %title(sprintf('Inner sum: %.2f, Pix: %.0f, Outer sum: %.2f, Pix: %.0f, Diff: %.2f',I2(ifile),masksize,I1(ifile),outersize,diffghostcheck(ifile)));
%     set(gca,'YDir','normal');
%     hold off;
%     if( (mod(ifile-1,subplot_tot) == (subplot_tot-1)) || (ifile == size(datafiles,1)) )
%         subplot_cntr = subplot_cntr + 1; % Increment
%         ext = '.png';
% %         imagename = [paths.ghostdiffcompdir,'set',num2str(subplot_cntr),'_',num2str(ifile-subplot_tot+1),'-',num2str(ifile),'_',num2str(subplot_setup(1)),'x',num2str(subplot_setup(2)),ext];
%         imagename = [paths.ghostdiffcompdir,'lim',num2str(min(valu_lims)),'to',num2str(max(valu_lims)),'_set',num2str(subplot_cntr),'_',num2str(ifile-subplot_tot+1),'-',num2str(ifile),'_',num2str(subplot_setup(1)),'x',num2str(subplot_setup(2)),ext];
%         print(h,imagename, '-dpng');
%     end
    
    %----- Single Image Plot Saver -----
%     h = figure(1);
%     clf;
%     set(h,'visible','off');
%     imagesc(data.image.calimage.*~data.mask.onemask);
%     a = colorbar;
%     a.Label.String = 'Intensity [nW m^{-2} sr^{-1}]';
%     caxis([-250,250]);
%     hold on;
%     pbaspect([1 1 1]);
%     xlabel('LORRI X Pixels');
%     ylabel('LORRI Y Pixels');
% %     th = 0:pi/50:2*pi;
% %     xunit = 58.2*cos(th)+128;
% %     yunit = 58.2*sin(th)+128;
% %     plot(xunit,yunit,'r');
% %     title(sprintf('Inner sum: %.2f, Pix: %.0f, Outer sum: %.2f, Pix: %.0f, Diff: %.2f',I2(ifile),masksize,I1(ifile),outersize,diffghostcheck(ifile)));
%     title(sprintf('Mean: %.2f Max: %.2f Min: %.2f',im_hist(ifile),max(im),min(im)));
%     set(gca,'YDir','normal');
%     ext = '.png';
% %     imagename = sprintf('%s%s%s',paths.ghostdiffcompdir,data.header.timestamp,ext);
%     imagename = sprintf('%s%s%s',paths.maskdir,data.header.timestamp,ext);
%     print(h,imagename, '-dpng');

    %% 
    %print field number
    field(ifile) = data.header.fieldnum;
    disp(data.header.fieldnum)
    date(ifile) = data.header.date_jd-data.header.launch_jd; %days from launch
    im_num(ifile) = ifile; %image number in sequence
    ccdtemp(ifile) = data.astrom.ccdtemp;
    fpubtemp(ifile) = data.astrom.fpubtemp;
    date_jd(ifile) = data.header.date_jd;
    bad(ifile) = data.header.bad;
    
    % Stop on a certain file
%     if ifile == 65 | ifile == 139 | ifile == 252 | ifile == 281 | ifile == 341
%         fprintf('ahhhh')
%         data.header.bad = 1;
%     else
%         data.header.bad = 0;
%     end
    
    % Stop on a certain field
%     if data.header.fieldnum == 7
%         fprintf('ahhhh')
%     end
    
end

%%
% th = 0:pi/50:2*pi;
% xunit = 58.2*cos(th)+128;
% yunit = 58.2*sin(th)+128;
% plot(xunit,yunit);
% set(gca,'YDir','normal');

figcnt = 1;

% Plot CCD temp vs. julian date to see if bad images are first in sequence
% figure(figcnt)
% colormap(jet)
% scatter(date_jd(field>1 & bad<1),ccdtemp(field>1 & bad < 1),[],field(field>1 & bad < 1))
% hold on
% scatter(date_jd(field>1 & bad>0),ccdtemp(field>1 & bad > 0),300,field(field>1 & bad > 0),'x')
% g = colorbar;
% ylabel(g,'Field number')
% xlabel('Julian Date')
% ylabel('CCD Temp [K]')
% figcnt = figcnt + 1;

% Plot FPUB temp vs. julian date to see if bad images are first in sequence
% figure(figcnt)
% colormap(jet)
% scatter(date_jd(field>1),fpubtemp(field>1),[],field(field>1))
% hold on
% scatter(date_jd(field>1 & bad>0),fpubtemp(field>1 & bad > 0),300,field(field>1 & bad > 0),'x')
% g = colorbar;
% ylabel(g,'Field number')
% xlabel('Julian Date')
% ylabel('FPUB Temp [K]')
% figcnt = figcnt + 1;

% Plot difference in mask mean before and after jail bar correction
figure(figcnt)
jailbarPlotter =  zeros(size(datafiles));
jailbarPlotter(jailbar_type > 1) = 1;
s1 = scatter(im_num(field>1 & bad < 1),(maskmean(field>1 & bad < 1)-maskmeanold(field>1 & bad < 1)).*exposure_time(field>1 & bad < 1),[],jailbarPlotter(field>1 & bad < 1));
hold on;
p1 = plot(im_num(field>1 & bad < 1),repmat(mean((maskmean(field>1 & bad < 1 & jailbarPlotter == 0)-maskmeanold(field>1 & bad < 1 & jailbarPlotter == 0)).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 0)),sum(field>1 & bad < 1),1),'--');
t1 = text(mean(im_num(field>1 & bad < 1)),0.20,num2str(mean((maskmean(field>1 & bad < 1 & jailbarPlotter == 0)-maskmeanold(field>1 & bad < 1 & jailbarPlotter == 0)).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 0))));
p2 = plot(im_num(field>1 & bad < 1),repmat(mean((maskmean(field>1 & bad < 1 & jailbarPlotter == 1)-maskmeanold(field>1 & bad < 1 & jailbarPlotter == 1)).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 1)),sum(field>1 & bad < 1),1),'--');
t2 = text(mean(im_num(field>1 & bad < 1)),-0.20,num2str(mean((maskmean(field>1 & bad < 1 & jailbarPlotter == 1)-maskmeanold(field>1 & bad < 1 & jailbarPlotter == 1)).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 1))));
xlabel('Image Number')
ylabel('Image Masked Mean Delta (JailBarCorrection - PreJailBar) [DN]')
colormap(jet);
colorbar;
figcnt = figcnt + 1;

% Plot mask mean (DN) before and after jail bar correction over all images
figure(figcnt)
s1_one = scatter(im_num(field>1 & bad < 1 & jailbarPlotter == 1)...
    ,maskmeanold(field>1 & bad < 1 & jailbarPlotter == 1).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 1)...
    ,[],'MarkerEdgeColor',[1 51/255 51/255]);
hold on;
s1_zero = scatter(im_num(field>1 & bad < 1 & jailbarPlotter == 0)...
    ,maskmeanold(field>1 & bad < 1 & jailbarPlotter == 0).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 0)...
    ,[],'MarkerEdgeColor',[0 102/255 204/255]);
% s2 = scatter(im_num(field>1 & bad < 1),maskmean(field>1 & bad < 1).*exposure_time(field>1 & bad < 1),'g','x');
xlabel('Image Number')
ylabel('Image Masked Mean [DN]')
legend([s1_one,s1_zero],{'Before Jail Bar Correction (High)','Before Jail Bar Correction (Low)'})
ylimz = ylim;
figcnt = figcnt + 1;

figure(figcnt)
s2_one = scatter(im_num(field>1 & bad < 1 & jailbarPlotter == 1)...
    ,maskmean(field>1 & bad < 1 & jailbarPlotter == 1).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 1)...
    ,[],'MarkerEdgeColor',[1 51/255 51/255]);
hold on;
s2_zero = scatter(im_num(field>1 & bad < 1 & jailbarPlotter == 0)...
    ,maskmean(field>1 & bad < 1 & jailbarPlotter == 0).*exposure_time(field>1 & bad < 1 & jailbarPlotter == 0)...
    ,[],'MarkerEdgeColor',[0 102/255 204/255]);
xlabel('Image Number')
ylabel('Image Masked Mean [DN]')
legend([s2_one,s2_zero],{'After Jail Bar Correction (High)','After Jail Bar Correction (Low)'})
ylim(ylimz); %make y axis match
figcnt = figcnt + 1;

% Masked image mean before and after jail bar correction colored by field
% number
figure(figcnt)
s1 = scatter(im_num(field>1 & bad < 1),maskmeanold(field>1 & bad < 1));
% hold on;
s2 = scatter(im_num(field>1 & bad < 1),maskmean(field>1 & bad < 1).*exposure_time(field>1 & bad < 1),[],field(field>1 & bad < 1));
xlabel('Image Number')
ylabel('Image Masked Mean [DN]')
legend([s1,s1],{'Before Jail Bar Correction','After Jail Bar Correction'})
figcnt = figcnt + 1;

% Plot comparison of diff ghost sub from model and data
figure(figcnt)
colormap(lines(length(min(field(field>1 & bad < 1)):1:max(field(field>1 & bad < 1)))))
scatter(diffghostcheck(field>1 & bad < 1),diffghostreal(field>1 & bad < 1),[],field(field>1 & bad < 1))
hold on;
% scatter(diffghostcheck(field>1 & bad > 0),diffghostreal(field>1 & bad > 0),300,field(field>1 & bad > 0),'x')
h = colorbar('Ticks',min(field(field>1 & bad < 1)):1:max(field(field>1 & bad < 1)));
ylabel(h, 'Field number')
xlabel('\lambdaI_{\lambda}^{G,D} [nW m^{-2} sr^{-1}]')
ylabel('\lambdaI_{\lambda}^{G,M} [nW m^{-2} sr^{-1}]')
[fitobject,gof,output] = fit(diffghostcheck,diffghostreal,'poly1');
xfit=linspace(min(diffghostcheck),max(diffghostcheck));
yfit=(fitobject.p1*xfit + fitobject.p2);
fit2 = plot(xfit,yfit);
title(sprintf('Fit: y = %.3fx + %.3f',fitobject.p1,fitobject.p2));
figcnt = figcnt + 1;

figure(figcnt)
% Plot hist of x-axis (diff ghost check from data) to see distribution
histogram(diffghostcheck(field>1),25)
xlabel('\lambdaI_{\lambda}^{G,D} [nW m^{-2} sr^{-1}]')
ylabel('N')
% Mean and std give acceptable range for y-int to be within 0
mean(diffghostcheck(field>1))
std(diffghostcheck(field>1))
figcnt = figcnt + 1;

% Plot comparison of hist peak of masked calibrated image and ref bias
figure(figcnt)
colormap(jet)
% Uncorrected most prob value and ref bias
scatter(im_ref(field>1 & bad < 1),im_hist(field>1 & bad < 1),[],field(field>1 & bad < 1))
% Corrected most prob value and ref bias
% scatter(im_ref(field>1),im_mostprob(field>1),[],field(field>1))
% Corrected mean and ref bias
% scatter(im_ref(field>1 & bad < 1),im_mean(field>1 & bad < 1),[],field(field>1 & bad < 1))
% hold on;
% scatter(im_ref(field>1 & bad > 0),im_mean(field>1 & bad > 0),300,field(field>1 & bad > 0),'x')
g = colorbar;
ylabel(g,'Field number')
xlabel('Reference Pixel Bias [DN]')
ylabel('Masked Image Mean [nW m^{-2} sr^{-1}]')
% ylabel('Corrected Most Probable Masked Image Value [nW m^{-2} sr^{-1}]')
% ylabel('Corrected Masked Image Mean [nW m^{-2} sr^{-1}]')