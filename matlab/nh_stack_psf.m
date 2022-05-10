function nh_stack_psf()
use_gaia = 1;

resize = 10;
pdim = 101;

psfx = ceil(pdim/2)-1;
psfy = ceil(pdim/2)-1;

% paths = get_paths_new();
% paths = get_paths_old();
paths = get_paths_lauer();

datafiles = dir(sprintf('%s*.mat',paths.datadir));

nfiles = size(datafiles,1);

for ifile=1:nfiles
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
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
    
    map = data.data;
    
    map = imresize(map,resize,'nearest');
%     map = RegridderZen(map,size(map)*resize); %better, flux conserving
    
    stack = zeros(pdim);
    shits = zeros(pdim);
    halfstack1 = zeros(pdim);
    halfshits1 = zeros(pdim);
    halfstack2 = zeros(pdim);
    halfshits2 = zeros(pdim);
    
    mymags = NaN.*ones(ncat,1);
    
    % loop over each catalog entry;
    for row = 1:ncat
        
        if use_gaia == 0
            
            mags = [B1mag(row),B2mag(row),R1mag(row),R2mag(row),I2mag(row)];
            
            sg = [B1sg(row),B2sg(row),R1sg(row),R2sg(row),I2sg(row)];
            
            lambda_mag = [425,462.5,645,650,810];
            
            whpl = (mags < 20) & (mags > 1);
            
            if sum(whpl) > 1
                thismag = nh_synthetic_photometry(lambda_mag(whpl),mags(whpl),'V');
            else
                thismag = NaN;
            end
            
        elseif use_gaia == 1
            
            thismag = Gmag(row);
        end
        
        %find x/y coordinate of the object
        [ypix, xpix] = radec2pix(RA(row),DEC(row), data.astrom);
        
        pixx = round(resize*xpix) - 9;
        pixy = round(resize*ypix) - 9;
        
        if ~isnan(thismag) & thismag < 16 & ...
                round(pixx-psfx) >= 1 & round(pixx+psfx) <= 2560 & ...
                round(pixy-psfy) >= 1 & round(pixy+psfy) <= 2560
            
            mymags(row) = thismag;
            
            thisx = [round(pixx-psfx):round(pixx+psfx)];
            thisy = [round(pixy-psfy):round(pixy+psfy)];
            
            if max(map(thisy,thisx)) < 2000
                stack = stack + map(thisy,thisx);
                
                shits = shits + ones(pdim);
                
                if randi(2) == 2
                    halfstack1 = halfstack1 + map(thisy,thisx);
                    halfshits1 = halfshits1 + ones(pdim);
                else
                    halfstack2 = halfstack2 + map(thisy,thisx);
                    halfshits2 = halfshits2 + ones(pdim);
                end
                
                %figure(1);
                %whpl = shits > 0;
                %tstack = zeros(pdim);
                %tstack(whpl) = stack(whpl) ./ shits(whpl);
                %imagesc(tstack)
                %colorbar
                %drawnow
                %pause(0.1)
            end
        end
        
    end
    
    stack = stack ./ shits;
    
    halfstack1 = halfstack1 ./ halfshits1;
    halfstack2 = halfstack2 ./ halfshits2;
    
    %% fit a gaussian to the grid
    [x y] = meshgrid(1:pdim,1:pdim);
    ft = fit_2D_gaussian(x,y,stack);
    ft1 = fit_2D_gaussian(x,y,halfstack1);
    ft2 = fit_2D_gaussian(x,y,halfstack2);
    
    data2 = load('/data/symons/NH_old_data/mat/2453978.8117112.mat');
    
    zoommap = data2.data.data(17:41,23:47);
    zoommap(7:8,3:4) = median(zoommap(:));
    
    zoomsize = size(zoommap);
    
    zoommap(round(zoomsize(1)./2)-1:round(zoomsize(1)./2)+1,...
        round(zoomsize(2)./2)+1:end)=...
        fliplr(zoommap(round(zoomsize(1)./2)-1:round(zoomsize(1)./2)+1,...
        1:round(zoomsize(2)./2)-1));
    
    zoommap = imresize(zoommap,resize,'nearest');
%     zoommap = RegridderZen(zoommap,size(zoommap)*resize);
    
    %figure(1); clf
    %plot(0.03.*zoommap(125,:)./max(zoommap(:)))
    %hold on
    %plot([1:numel(stack(58,:))]+69,stack(58,:)./max(stack(:)))
    
    [xgrid,ygrid] = meshgrid([-psfx:psfx]);
    
    rad = sqrt((xgrid-5).^2 + (ygrid-5).^2);
    
    mask = zeros(pdim);
    
    whpl = rad <= 19;
    mask(whpl) = 1;
    whouter = rad > 40;% & rad < 21;
    stack = stack - mean(stack(whouter));
    
    whpl = rad <= 10;
    
    newpsf = 0.05.*zoommap./max(zoommap(:));
    newstack = stack ./ max(stack(:));
    
    [xgrid1,ygrid1] = meshgrid([-psfx:psfx]);
    rad1 = sqrt((xgrid1-5).^2 + (ygrid1-5).^2);
    
    whpl1 = rad1 < 15.5;
    
    [xgrid2,ygrid2] = meshgrid([1:size(zoommap,1)]);
    rad2 = sqrt((xgrid2-125).^2 + (ygrid2-127).^2);
    
    whpl2 = rad2 < 15.5;
    
    newpsf(whpl2) = newstack(whpl1);
    
    %figure(1); clf;
    %imagesc(newpsf)
    %caxis([0,0.1])
    %colorbar
    
    tempint = newpsf./max(newpsf(:));
    disp(sprintf('PSF Omega Factor is: %5.3f',...
        sum(tempint(:)) ./ (pdim./resize).^2));
    
    newpsfp = wshift('2D',newpsf,[3,0]);
    tenxpsf = imresize(wshift('2D',newpsf,[3,1]),1./(pdim./(size(zoommap,1)./resize)),'nearest');
%     tenxpsf = RegridderZen(wshift('2D',newpsf,[3,1]),ceil(size(newpsf)*1./(pdim./(size(zoommap,1)./resize))));
    tenxpsf = tenxpsf ./ max(tenxpsf(:));
    
    fourxpsf = imresize(wshift('2D',newpsf,[1,-1]),1./resize,'nearest');
%     fourxpsf = RegridderZen(wshift('2D',newpsf,[1,-1]),ceil(size(newpsf)*1./resize));
    fourxpsf = fourxpsf ./ max(fourxpsf(:));
    
    thispsf = imresize(wshift('2D',newpsf,[-6,-15]),1./(pdim./(size(zoommap,1)./resize).*resize),'nearest');
%     thispsf = RegridderZen(wshift('2D',newpsf,[-6,-15]),ceil(size(newpsf)*1./(pdim./(size(zoommap,1)./resize).*resize)));
    
    thispsfp = thispsf ./ sum(thispsf(:));
    
%     figure(2); clf;
%     imagesc(thispsf)
%     colorbar
    
    modelpsf = tenxpsf(1:61,1:61);
    modelpsf = modelpsf ./ sum(modelpsf(:));
    
    psf.psf = thispsfp;
    psf.modelpsf = modelpsf;
    psf.onexpsf = thispsf ./ max(thispsf(:));
    psf.fourxpsf = fourxpsf;
    psf.tenxpsf = tenxpsf;
    psf.fourtyxpsf = newpsfp./max(newpsfp(:));
    psf.pixelwidths = [4.3,0.43,4.3,4.3./4.04,0.43,4.3./40.4];
    psf.centerpix = [4,31,4,13,31,125];
    
    data.psf = psf;
    
    save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
    
    %save('lookup/nh_lorri_psf.mat','psf');
    
end
end

