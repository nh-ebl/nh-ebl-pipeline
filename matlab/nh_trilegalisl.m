function data = nh_trilegalisl(paths, data, tri_gaia, tri_mag, max_mag, save_file)

fieldnum = data.header.fieldnum;

% step 1: compute the area of a LORRI image and pull out that many stars
% from the list at random
surveyarea = data.astrom.imagew.*data.astrom.imageh.*...
    data.cal.pixsize_arcsec.^2 ./ 3600.^2;

% step 2: read in the trilegal information
load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/isltrilegal_%02d.mat',fieldnum));

nfiles = numel(V);

isltot = zeros(nfiles,1);
islmasked = zeros(nfiles,1);

for jfile=1:nfiles
    
    mag = V(jfile).mlf;
    
    % step 3: convert from LORRI-band mag to flux
    Fcat = data.cal.vzero .* 10.^(-V(jfile).mlf/2.5);
    
    % step 4: make a mask function
    if tri_gaia == 0
        whpl = V(jfile).V > max_mag;
    elseif tri_gaia == 1
        whpl = V(jfile).V > tri_mag;
    end
    
    magcat = mag(whpl);
    
    % step 5: convert to surface brightness
    lIltot = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea .* (pi./180).^2);
    lIlcat = 1e-26.*1e9.*data.cal.nu.*Fcat(whpl)./(surveyarea .* (pi./180).^2);
    
    isltot(jfile) = sum(lIltot);
    islmasked(jfile) = sum(lIlcat);
    
end

star_list = horzcat(magcat,(lIlcat.*(data.astrom.imagew.*data.astrom.imageh)));

if tri_gaia == 0
    filename = strcat(num2str(max_mag),'_tri_only.mat');
elseif tri_gaia == 1
    filename = strcat(num2str(tri_mag),'_tri.mat');
end

if save_file == 1
    save(filename,'star_list');
end

islout.isltotmean = mean(isltot);
islout.isltoterr = std(isltot);

islout.islmaskedmean = mean(islmasked);
islout.islmaskederr = std(islmasked);

% fileout = sprintf('%sisl/%s',paths.tridir,datafiles(ifile).name);
%
% save(fileout,'islout')

data.isl.tritotmean = islout.isltotmean;
data.isl.tritoterr = islout.isltoterr;
data.isl.trimean = islout.islmaskedmean;
data.isl.trierr = islout.islmaskederr;

end


