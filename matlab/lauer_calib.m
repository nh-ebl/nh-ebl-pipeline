clear variables
close all

%get paths for new data files or old data files
paths = get_paths_lauer();
datadir = paths.datadir;

%desired image and location
datafile = '2458731.2928864.mat';
sourcex = 73.3068 + 1;
sourcey = 58.6167 + 1;

%load individual mat file for image
load(sprintf('%s%s',datadir,datafile));
exptime = data.header.exptime;

%do aperture photometry
skysub_apsum = ap_photom_calib(data.image.rawimage.*data.header.exptime, data.image.rawimage.*~data.mask.onemask.*data.header.exptime, sourcex, sourcey, 3, 6, 9);

fprintf('done');