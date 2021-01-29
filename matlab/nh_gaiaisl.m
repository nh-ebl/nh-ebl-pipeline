function data = nh_gaiaisl(data,paths,tri_mag,max_mag,save_file)

fieldnum = data.header.fieldnum;

% step 1: compute the area of a LORRI image and pull out that many stars
% from the list at random
surveyarea = data.astrom.imagew.*data.astrom.imageh.*...
    data.cal.pixsize_arcsec.^2 ./ 3600.^2;

% step 2: load gaia catalog file
load(sprintf('%smat_files/field_%d_data.mat',paths.gaiadir,fieldnum));

mag = Gmag;

% step 3: convert from LORRI-band mag to flux
Fcat = data.cal.vzero .* 10.^(-Gmag/2.5);

% step 4: make a mask function
whpl = (Gmag > max_mag) & (Gmag < tri_mag);

magcat = mag(whpl);

% step 5: convert to surface brightness
lIltot = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea .* (pi./180).^2);
lIlcat = 1e-26.*1e9.*data.cal.nu.*Fcat(whpl)./(surveyarea .* (pi./180).^2);

star_list = horzcat(magcat,(lIlcat.*(data.astrom.imagew.*data.astrom.imageh)));

if save_file == 1
    filename = strcat(num2str(max_mag),'_gaia.mat');
    save(filename,'star_list');
end

isltot = sum(lIltot);
islmasked = sum(lIlcat);

data.isl.gaiatot = isltot;
% data.isl.tritoterr = triout.isltoterr;
data.isl.gaiamean = islmasked;
% data.isl.trierr = triout.islmaskederr;

end