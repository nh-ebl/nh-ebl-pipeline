%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_pipeline.m
%%  Jun 2016, MZ
%%  For each file in the NH data directories, this program:
%%   1) Reads it in.
%%   2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nh_pipeline()
clear all
close all
%get paths for new data files or old data files
% paths = get_paths_old_ghosts();
% paths = get_paths_old();
paths = get_paths_new();
datadir = paths.datadir;

%% Set user parameters here

%set portion of pipeline to run
procstepflag = 1; 

%use USNOB1 or Gaia for star masking
use_gaia = 1; %0 for USNOB1
new_star_mask = 0; %0 for skipping star masking and using existing mask

%use Gaia or just Trilegal for ISL calc
tri_gaia = 0; %0 for only Trilegal, 1 for Gaia and Trilegal

%keep or omit galaxies from Gaia data
gals = 0; %0 for galaxies removed, 1 for galaxies included

%set max mag for masking (mask all stars brighter than this), also
%calculate gaia isl between this and 20
max_mag = 21 ; %17.75 old value

%set min mag for trilegal isl (calculate isl for stars fainter than this)
tri_mag = 20;

%set max mag for psf wing isl calc (from USNOB1 or Gaia based on masking catalog)
wing_mag = max_mag;

%optionally save lists of stars and flux
save_file = 0; %1 to save

%%
%load directory of saved mat files for all images
datafiles = dir(sprintf('%s*.mat',datadir));

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


totalghosts = 0;
totalrealghosts = 0;

if procstepflag == 3
    evenMinusOdd = zeros( size(datafiles,1), 1); %preallocate
    evenMinusOddFixed = zeros( size(datafiles,1), 1); %preallocate
    oddMinusEvenref = zeros( size(datafiles,1), 1);
    oddMinusEvenrefFixed = zeros( size(datafiles,1), 1);
end

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

subplot_setup = [3,4]; % [# Rows, # Columns] for subplot
subplot_tot = prod(subplot_setup); % Number of plots per subplot
subplot_cntr = 0; % Counts for the subplot
valu_lims = [-250,250]; % Limits the plot dynamic range
origCmap = parula(4096); % Based on https://www.mathworks.com/matlabcentral/answers/307318-does-matlab-have-a-nonlinear-colormap-how-do-i-make-one#comment_695509
dataMax = max(valu_lims);
dataMin = min(valu_lims);
centerPoint = 0;
scalingIntensity = 6;
x = 1:length(origCmap); 
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*2047/max(x)+1; 
newCmap = interp1(x, origCmap, 1:2048);
newCmap(end,:) = [1,1,1]; % Make last value in cmap white (yellow end)
newCmap(1,:) = [1,1,1]; % Make 1st valuein cmap white (blue end)

%loop through all images
for ifile=1:size(datafiles,1)
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    %load individual mat file for image
    load(sprintf('%s%s',datadir,datafiles(ifile).name));
    disp(data.header.rawfile)
    
    % Plot new ghost mask - im>0 for all masks = 1, im for weighted by
    % repeated location
%     im = im + data.mask.ghostmask;
%     imagesc(im>0);
    
    % Calculate distance from ghost to center of fov + ghost radius
    ghostdist(ifile) = data.ghost.ghostdistcent + 21.5;
    
    % Calculate check of diffuse ghost sum from data (ghost region - not
    % ghost region)
    masksize = sum(data.mask.onemask.*ghostdiffmask,'all'); % Mask-based number of pixels
    outersize = sum(~data.mask.onemask.*~ghostdiffmask,'all'); % Mask-based number of pixels
    I2(ifile) = 1/masksize*sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all');
    I1(ifile) = 1/outersize*sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all');
    diffghostcheck(ifile) = 1/masksize*sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all')-1/outersize*sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all');
%     diffghostcheck(ifile) = 1/(256*256)*(1/masksize*sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all')-1/outersize*sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all'));
%     diffghostcheck(ifile) = 1/(256*256)*(sum(data.image.calimage.*~data.mask.onemask.*ghostdiffmask,'all')-sum(data.image.calimage.*~data.mask.onemask.*~ghostdiffmask,'all'));
    
    % Diffuse ghost contribution from model
    diffghostreal(ifile) = data.ghost.diffusesub;
%     diffghostreal(ifile) = data.ghost.diffusesub/masksize;

    % Plot masked, calibrated image with red circle indicating ghost region
    % Title includes sum and npix of inner and outer regions, and difference of
    % inner - outer, which estimates diff ghosts
    
    %----- Multi-Image Plot Saver -----
    if( mod(ifile-1,subplot_tot) == 0 )        
        h = figure(1);
        if( ifile == 1 )
            set(h, 'units', 'normalized'); % Force full screen
            set(h, 'outerposition', [0 0 1 1]); % Force full screen
%             set(h,'visible','off'); % Causes figure sizing problems
        end
        clf;
        % tiledlayout instead of subplot per https://www.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html
        tiledlayout(h,subplot_setup(1),subplot_setup(2),'TileSpacing','Compact','Padding','Compact'); 
    end
%     subplot(subplot_setup(1),subplot_setup(2),mod(ifile-1,subplot_tot)+1); % Tiled layout instead
    nexttile; % Tiled layout instead
    imagesc(data.image.calimage.*~data.mask.onemask);
    a = colorbar;
    a.Label.String = 'Intensity [nW m^{-2} sr^{-1}]';
    caxis([min(valu_lims),max(valu_lims)]);
%     colormap(newCmap);
    hold on;
    pbaspect([1 1 1]);
    xlabel('LORRI X Pixels');
    ylabel('LORRI Y Pixels');
    th = 0:pi/50:2*pi;
    xunit = 58.2*cos(th)+128;
    yunit = 58.2*sin(th)+128;
    plot(xunit,yunit,'r');
    title([num2str(ifile),' ',num2str(data.header.timestamp)])
    %title(sprintf('Inner sum: %.2f, Pix: %.0f, Outer sum: %.2f, Pix: %.0f, Diff: %.2f',I2(ifile),masksize,I1(ifile),outersize,diffghostcheck(ifile)));
    set(gca,'YDir','normal');
    hold off;
    if( (mod(ifile-1,subplot_tot) == (subplot_tot-1)) || (ifile == size(datafiles,1)) )
        subplot_cntr = subplot_cntr + 1; % Increment
        ext = '.png';
%         imagename = [paths.ghostdiffcompdir,'set',num2str(subplot_cntr),'_',num2str(ifile-subplot_tot+1),'-',num2str(ifile),'_',num2str(subplot_setup(1)),'x',num2str(subplot_setup(2)),ext];
        imagename = [paths.ghostdiffcompdir,'lim',num2str(min(valu_lims)),'to',num2str(max(valu_lims)),'_set',num2str(subplot_cntr),'_',num2str(ifile-subplot_tot+1),'-',num2str(ifile),'_',num2str(subplot_setup(1)),'x',num2str(subplot_setup(2)),ext];
        print(h,imagename, '-dpng');
    end
    
    %----- Single Image Plot Saver -----
%     h = figure(1);
%     clf;
%     set(h,'visible','off');
%     imagesc(data.image.calimage.*~data.mask.onemask);
%     a = colorbar;
%     a.Label.String = 'Intensity [nW m^{-2} sr^{-1}]';
% %     caxis([-5,5]);
%     hold on;
%     pbaspect([1 1 1]);
%     xlabel('LORRI X Pixels');
%     ylabel('LORRI Y Pixels');
%     th = 0:pi/50:2*pi;
%     xunit = 58.2*cos(th)+128;
%     yunit = 58.2*sin(th)+128;
%     plot(xunit,yunit,'r');
%     title(sprintf('Inner sum: %.2f, Pix: %.0f, Outer sum: %.2f, Pix: %.0f, Diff: %.2f',I2(ifile),masksize,I1(ifile),outersize,diffghostcheck(ifile)));
%     set(gca,'YDir','normal');
%     ext = '.png';
%     imagename = sprintf('%s%s%s',paths.ghostdiffcompdir,data.header.timestamp,ext);
%     print(h,imagename, '-dpng');

    %print field number
    disp(data.header.fieldnum)
    
    if ifile == 55
        fprintf('ahhhh')
    end
    
    if procstepflag == 1
        fprintf('Masking and ghost analysis');
        %% Masking and optical ghosts
        %save ra, dec, pix, and mag of stars in image that could be causing ghosts
        %         [data, ghostcount, realghostcount] = nh_findghoststar(data,paths,use_gaia);
        %                 totalghosts = totalghosts + ghostcount;
        %                 totalrealghosts = totalrealghosts + realghostcount;
        
        %
        %manually mask portions of image
        %                 data = nh_make_manmask(data,paths);
        
        %mask out stars and ghosts in image (also stat and clip mask)
        %may need to redo catalog depending on gals
        %                 catalog_data_gaia(gals,paths)
%                         data = nh_makemask(data,paths,3,use_gaia,new_star_mask, max_mag, save_file);
        
        %         maskim = data.data.*~data.mask.onemask;
        %         maskim(maskim==0) = NaN;
        %         checkmax(ifile) = nanmax(nanmax(maskim));
        %         checkmin(ifile) = nanmin(nanmin(maskim));
        
        %print total masked surface brightness (sum of all masked pixels)
        %         data.cal.sbconv .* sum(sum(data.data(data.mask.onemask)./ data.astrom.exptime))
        
        %print mask fraction
        %         data.mask.maskfrac
        
        %overwrite data file with changes
%                 save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 2
        fprintf('Meta data and dark current');
        %% Meta data and dark current
        %write header values, constants, coordinates, and dark reference pixel data to data file
        %(uses mask, may need to redo mask first)
        %changes values used in calibration
        data = nh_add_meta(data);
        
        %overwrite data file with changes
        save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 3
        fprintf('Jail bars and calibration');
        %% Jail bar correction and calibration
        %perform jail bar correction on data.data and data.ref - need mask first, but can redo from original data
        %                 [data, evenMinusOdd(ifile), evenMinusOddFixed(ifile), oddMinusEvenref(ifile), oddMinusEvenrefFixed(ifile)] = nh_jail_bars(data,paths);
        
        %save calibrated image in surface brightness units (nh_calcref saves refcorr which is used here)
        %must redo calibration after jail bar correction, but nh_calcref must be run first
        data = nh_calibrate(data,paths);
        
        %overwrite data file with changes
        save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 4
        fprintf('ISL calculation');
        %% ISL calculation
        
        %calculate ISL from USNOB1 and Trilegal
        %may need to redo catalog depending on gals
        %         catalog_data_gaia(gals,paths)
        data = nh_calcisl(data, paths, use_gaia, tri_gaia, tri_mag, wing_mag, max_mag, save_file);
        
        %         wing = data.isl.usnowing
        %         triisl = data.isl.trimeanmasksize
        %         gaiaisl = data.isl.gaiamean
        
        %overwrite data file with changes
        %         save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');
        
    end
    
    if procstepflag == 5
        fprintf('DGL calculation');
        %% DGL calculation
        %         if data.header.fieldnum == 8
        %             fprintf('its here');
        %         end
        
        %calculate DGL using Planck (or IRIS) data
        data = nh_calcdgl(data, paths);
        
        %overwrite data file with changes
        save(sprintf('%s%s',datadir,datafiles(ifile).name),'data');
        
    end
    
    %     mydate(ifile) = data.header.date_jd;
    %     mytemp(ifile) = data.header.ccdtemp;
    %     mymean(ifile) = data.stats.calmean;
    %     myref(ifile) = data.ref.mean;
    %     myeng(ifile) = data.ref.engmean;
    %     %myisl(ifile) = data.isl.trimean;
    %
    %     if 0
    %     figure(1); clf
    %     imagesc(data.image.calimage)
    %     caxis([0,3000])
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	'(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_raw.png',data.header.timestamp));
    %     figure(2); clf
    %     imagesc(data.mask.mask)
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	  '(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_mask.png',data.header.timestamp));
    %     figure(3); clf
    %     imagesc(data.image.calimage.*~data.mask.mask)
    %     caxis([0,1500])
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	'(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_masked.png',data.header.timestamp));
    %       end
    
    
    %dbstop
    
    
end

% th = 0:pi/50:2*pi;
% xunit = 58.2*cos(th)+128;
% yunit = 58.2*sin(th)+128;
% plot(xunit,yunit);
% set(gca,'YDir','normal');

% scatter(diffghostcheck,diffghostreal)
% hold on;
% xlabel('\lambdaI_{\lambda}^{G,D} [nW m^{-2} sr^{-1}]')
% ylabel('\lambdaI_{\lambda}^{G,M} [nW m^{-2} sr^{-1}]')
% [fitobject,gof,output] = fit(diffghostcheck,diffghostreal,'poly1');
% xfit=linspace(min(diffghostcheck),max(diffghostcheck));
% yfit=(fitobject.p1*xfit + fitobject.p2);
% fit2 = plot(xfit,yfit);
% title(sprintf('Fit: y = %.3fx + %.3f',fitobject.p1,fitobject.p2));




totalghosts
totalrealghosts
max(ghostdist)

% plot(oddMinusEvenref);
% plot(oddMinusEvenrefFixed);
% plot(evenMinusOdd);
% plot(evenMinusOddFixed);



%figure(4);
%plot(mydate - data.header.launch_jd,mytemp,'o');


end
