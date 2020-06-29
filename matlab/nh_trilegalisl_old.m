function nh_trilegalisl()

paths = get_paths();

datafiles = dir(sprintf('%s*.mat',paths.datadir));

%warning('off','all')

for ifile=1:size(datafiles,1)
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    fieldnum = data.header.fieldnum;
    
    % step 1: compute the area of a LORRI image and pull out that many stars
    % from the list at random
    surveyarea = data.astrom.imagew.*data.astrom.imageh.*...
        data.cal.pixsize_arcsec.^2 ./ 3600.^2;
    
    if fieldnum == 5
        
        % step 2: read in the trilegal information
        load(sprintf('/home/symons/isl_trilegal/isltrilegal_%02d.mat',fieldnum));
        
        nfiles = numel(V);
        
        isltot = zeros(nfiles,1);
        islmasked = zeros(nfiles,1);
        
        for jfile=1:nfiles
            
            mag = V(jfile).mlf;
            
            % step 3: convert from LORRI-band mag to flux
            Fcat = 3055 .* 10.^(-V(jfile).mlf/2.5);
            
            % step 4: make a mask function
            whpl = V(jfile).V > data.mask.maxmag;
            
            magcat = mag(whpl);
            
            % step 5: convert to surface brightness
            lIltot = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea .* (pi./180).^2);
            lIlcat = 1e-26.*1e9.*data.cal.nu.*Fcat(whpl)./(surveyarea .* (pi./180).^2);
            
            isltot(jfile) = sum(lIltot);
            islmasked(jfile) = sum(lIlcat);
            
        end
        
        islout.isltotmean = mean(isltot);
        islout.isltoterr = std(isltot);
        
        islout.islmaskedmean = mean(islmasked);
        islout.islmaskederr = std(islmasked);
        
        fileout = sprintf('%sisl/%s',paths.tridir,datafiles(ifile).name);
        
        save(fileout,'islout')
        
    end
    
end

%warning('on','all')

% end
