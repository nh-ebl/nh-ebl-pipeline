function usnoisl = nh_usnoisl(data, paths, use_gaia, wing_mag, save_file, flag_method, errflag_psf)

load('run_params.mat','params')
if params.err_on == 1
    datastruct = data.(params.err_str);
    maskdir = data.err_mags;
else
    datastruct = data;
    maskdir = data;
end

resize = 10;

% If new method, load per data file psf
if strcmp(flag_method,'new') == 1
    thispsf = data.psf.modelpsf;
    % If old method, load single saved old psf
elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'old') == 1)
    load('lookup/nh_lorri_psf.mat');
    thispsf = psf.modelpsf;
end

starimage = zeros(256.*resize);
wingimage = zeros(256.*resize);

% load up the corresponding catalog file
if use_gaia == 1
    load(sprintf('%smat_files/field_%d_data.mat',paths.gaiadir,data.header.fieldnum));
    
    % figure out the length of the catalog
    [ncat,~] = size(RA);
elseif use_gaia == 0
    load(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum));
    
    % figure out the length of the catalog
    [~,ncat] = size(RA);
end

mymag = zeros(ncat,1);
myflux = zeros(ncat,1);
ypixlist = zeros(ncat,1);
xpixlist = zeros(ncat,1);

if ~exist('xpix','var') || ~exist('ypix','var')
    %only run if need to make xpix/ypix (should not run since makemask made them)
    ypix = zeros(ncat,1);
    xpix = zeros(ncat,1);
    % convert radec2pix parallel now since will be required (later loop
    % can't be paralleled w/o work)
    parpoolobj = gcp('nocreate'); % check for thread pool
    if isempty(parpoolobj)
        maxNumCompThreads(7); % Declare number of threads to use
        parpool('threads');
    else
        if ~isa(parpoolobj,'parallel.ThreadPool')
            delete(parpoolobj); %want threads here
            maxNumCompThreads(7); % Declare number of threads to use
            parpool('threads');
        end
    end
    RApar = RA;
    DECpar = DEC; %parfor was mad about these not ever being declared officially
    parfor row = 1:ncat
        %parallel for speed, later loop needs work to parallize
        [ypix(row), xpix(row)] = radec2pix(RApar(row),DECpar(row), data.astrom);
    end
    % if params.err_gals == 0
        % save the calc'd ypix/xpix the corresponding catalog file
        if use_gaia == 1
            save(sprintf('%smat_files/field_%d_data.mat',paths.gaiadir,data.header.fieldnum),'xpix','ypix','-append');
        elseif use_gaia == 0
            save(sprintf('%sfield_%d_data.mat',paths.catdir,data.header.fieldnum),'xpix','ypix','-append');
        end
    % elseif params.err_gals == 1
    %     % save the calc'd ypix/xpix the corresponding catalog file
    %     if use_gaia == 1
    %         save(sprintf('%smat_files/field_%d_mc/%i', paths.gaiadir,data.header.fieldnum,mc),'xpix','ypix','-append');
    %     end
    % end
end

% loop over each catalog entry;
for row = 1:ncat
    
    if use_gaia == 0
        
        mags = [B1mag(row),B2mag(row),R1mag(row),R2mag(row),I2mag(row)];
        sg = [B1sg(row),B2sg(row),R1sg(row),R2sg(row),I2sg(row)];
        lambda_mag = [425,462.5,645,650,810];
        whpl = (mags < 21) & (mags > 1);
        if sum(whpl) > 1
            thismag = nh_synthetic_photometry(lambda_mag(whpl),mags(whpl),'LORRI');
        else
            thismag = NaN;
        end
        
    elseif use_gaia == 1
        
        if errflag_psf == 1
            thismag = Gmag(row) + Gmagerr(row)*randn(1); %If including mag error, include up to the full reported Gaia error with Gaussian probability
        elseif errflag_psf == 0
            thismag = Gmag(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value
        end        
    end
    
    %find x/y coordinate of the object
    % [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom); %par now
    
    if ypix(row) >= 1 & ypix(row) <= 256 & xpix(row) >= 1 & xpix(row) <= 256 & ~isnan(thismag)
        
        thisflux = data.cal.vzero .* 10^(-thismag./2.5);
        
        starimage(round(ypix(row).*resize),round(xpix(row).*resize)) = ...
            starimage(round(ypix(row).*resize),round(xpix(row).*resize)) + thisflux;
        
        if thismag < wing_mag
            wingimage(round(ypix(row).*resize),round(xpix(row).*resize)) = ...
                wingimage(round(ypix(row).*resize),round(xpix(row).*resize)) + thisflux;
            mymag(row) = thismag;
            myflux(row) = sum(sum(wingimage(round(ypix(row).*resize),round(xpix(row).*resize))));
            xpixlist(row) = xpix(row);
            ypixlist(row) = ypix(row);
        end
        
    end
    
end

starimage_conv = conv2(starimage,thispsf,'same');
wingimage_conv = conv2(wingimage,thispsf,'same');

if strcmp(flag_method,'old') == 1
    starimage_dec = resize.*imresize(starimage_conv,1./resize,'nearest');
    wingimage_dec = resize.*imresize(wingimage_conv,1./resize,'nearest');
elseif (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
    starimage_dec = RegridderZen(starimage_conv,ceil(size(starimage_conv)*1./resize));
    wingimage_dec = RegridderZen(wingimage_conv,ceil(size(wingimage_conv)*1./resize));
end

starimage_cal = starimage_dec .* data.cal.nu .* 1e-26 .* 1e9 ./ ...
    data.cal.omega_pix;
wingimage_cal = wingimage_dec .* data.cal.nu .* 1e-26 .* 1e9 ./ ...
    data.cal.omega_pix;

wingimage_cal_masked = wingimage_cal;
wingimage_cal_masked(maskdir.mask.onemask) = nan;
psfsize = ceil(size(thispsf)/resize/2);

xpixlist = xpixlist(xpixlist~=0);
ypixlist = ypixlist(ypixlist~=0);
mycalflux = zeros(length(xpixlist),1);

for row = 1:length(xpixlist)
    if round(ypixlist(row))-psfsize(1) >= 1 & round(ypixlist(row))+psfsize(1) <= 256 & round(xpixlist(row))-psfsize(1) >= 1 & round(xpixlist(row))+psfsize(1) <= 256
        star = wingimage_cal_masked(round(ypixlist(row))-psfsize(1):round(ypixlist(row))+psfsize(1),round(xpixlist(row))-psfsize(1):round(xpixlist(row))+psfsize(1));
        mycalflux(row) = nansum(nansum(star));
    end
end

mymag = mymag(mymag~=0);
mymag = mymag(mycalflux~=0);
mycalflux = mycalflux(mycalflux~=0);

% myflux = myflux(myflux~=0);
% mycalflux = data.cal.nu .* 1e-26 .* 1e9 ./ data.cal.omega_pix .* myflux;

star_list = horzcat(mymag,mycalflux);

if save_file == 1
    filename = strcat(num2str(wing_mag),'_wing.mat');
    save(filename,'star_list');
end

usnoisl.isltot = mean(starimage_cal(~maskdir.mask.onemask));
usnoisl.totimage = starimage_cal;
usnoisl.islwing = mean(wingimage_cal(~maskdir.mask.onemask)); %this one gets used
usnoisl.wingimage = wingimage_cal;
usnoisl.islfaint = usnoisl.isltot - usnoisl.islwing;

% h = figure();
% clf;
% imagesc(wingimage_cal.*~maskdir.mask.onemask)
% set(h,'visible','off');
% a = colorbar;
% a.Label.String = 'Intensity [nW]';
% pbaspect([1 1 1]);
% xlabel('LORRI X Pixels');
% ylabel('LORRI Y Pixels');
% caxis([0,10]);
% title(sprintf('%8.3f',usnoisl.islwing));
% set(gca,'YDir','normal');
% ext = '.png';
% if not(isfolder([paths.wingdir]))
%     mkdir([paths.wingdir])
% end
% imagename = sprintf('%s%s%s',paths.wingdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

end