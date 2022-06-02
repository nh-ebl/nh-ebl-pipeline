function make_trilegalcats_gaia()

% paths = get_paths_new();
% paths = get_paths_old();
% paths = get_paths_lauer();
paths = get_paths_newest();

if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/good/') == 1
    old = 1;
    new = 0;
    lauer = 0;
    smallfields = [3,5,6,7]; % Old files
    nfields = 4;
elseif strcmp(paths.datadir,'/data/symons/nh_data/mat/') == 1
    new = 1;
    old = 0;
    lauer = 0;
    smallfields = [1,5,6,7,8]; % New files
    nfields = 5;
elseif strcmp(paths.datadir,'/data/symons/nh_data_lauer/mat/') == 1
    new = 0;
    old = 0;
    lauer = 1;
    smallfields = [1,2,3,4,5,6,7]; % New files
    nfields = 7;
elseif strcmp(paths.datadir,'/data/symons/nh_data_new/mat/') == 1
    new = 0;
    old = 0;
    lauer = 0;
    newest = 1;
    smallfields = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % New files
    nfields = 23;
end

%% read in trilegal dat files

for ifield=1:nfields %nfields % Make this be any set of numerically sequential fields - otherwise do one at a time

    disp(sprintf('On field %d of %d.',ifield,nfields));

    indir = dir(sprintf('%s/gaia/%02d/*.dat',paths.tridir,ifield));

    nfiles = numel(indir);

    for ifile=1:nfiles

        disp(sprintf('On file %d of %d.',ifile,nfiles));

        disp(sprintf('Reading file...'));

        % step 1: read in the trilegal catalog information
        filename = sprintf('%s/gaia/%02d/%s',paths.tridir,ifield,...
            indir(ifile).name);

        [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,G,a11,a12,a13] = ...
            textread(filename,['%d %f %f %f %f %f %f %f %f %f %f '...
            '%f %f %f %f'],'commentstyle','shell');

        clear a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13;

        nstars = numel(G);

        whpl = ifield == smallfields;
        if sum(whpl) == 0 % This skips currently, would need to be 1 to not skip
            nrands = nstars;
            nstars = round(0.0847 .* nstars); %should be 0.0841 to be 0.29*0.29 -> can convert 1 deg^2 of stars to 0.0841 deg^2
            ind = randperm(nrands,nstars);
            G = G(ind);
        end

        disp(sprintf('Calculating magnitudes...'));

        Vp(ifile).G = G;

    end

    V = Vp;

    if new == 1
        save(sprintf('lookup/trilegal/new_data/gaia/isltrilegal_%02d.mat',ifield),'V','-v7.3'); %Save location for new files
    elseif old == 1
        save(sprintf('/home/symons/isl_trilegal/gaia/isltrilegal_%02d.mat',ifield),'V','-v7.3'); % Save location for old files
    elseif lauer == 1
        save(sprintf('lookup/trilegal/lauer_data/gaia/isltrilegal_%02d.mat',ifield),'V','-v7.3'); % Save location for lauer files
    elseif newest == 1
        save(sprintf('lookup/trilegal/newest_data/gaia/isltrilegal_%02d.mat',ifield),'V','-v7.3'); % Save location for newest files
    end

end

end

