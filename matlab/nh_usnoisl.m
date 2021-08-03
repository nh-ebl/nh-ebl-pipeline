function usnoisl = nh_usnoisl(data, paths, use_gaia, wing_mag, save_file, flag_method)

resize = 10;

datafiles = dir(sprintf('%s*.mat',paths.datadir));

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
        
        thismag = Gmag(row); %+ randn(1) .* 0.25; %need to know what is possible gaia mag error to change this value
        
    end
    
    %find x/y coordinate of the object
    [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
    
    if ypix >= 1 & ypix <= 256 & xpix >= 1 & xpix <= 256 & ~isnan(thismag)
        
        thisflux = data.cal.vzero .* 10^(-thismag./2.5);
        
        starimage(round(ypix.*resize),round(xpix.*resize)) = ...
            starimage(round(ypix.*resize),round(xpix.*resize)) + thisflux;
        
        if thismag < wing_mag
            wingimage(round(ypix.*resize),round(xpix.*resize)) = ...
                wingimage(round(ypix.*resize),round(xpix.*resize)) + thisflux;
            mymag(row) = thismag;
            myflux(row) = sum(sum(wingimage(round(ypix.*resize),round(xpix.*resize))));
            xpixlist(row) = xpix;
            ypixlist(row) = ypix;
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
wingimage_cal_masked(data.mask.onemask) = nan;
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

usnoisl.isltot = mean(starimage_cal(~data.mask.onemask));
usnoisl.totimage = starimage_cal;
usnoisl.islwing = mean(wingimage_cal(~data.mask.onemask)); %this one gets used
usnoisl.wingimage = wingimage_cal;
usnoisl.islfaint = usnoisl.isltot - usnoisl.islwing;

%   figure(4); clf
%   imagesc(wingimage_cal.*~data.mask.onemask)
%   title(sprintf('%8.3f',usnoisl.islwing));
%   colorbar
%   caxis([0,10])
%   drawnow

end