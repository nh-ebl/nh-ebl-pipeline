function data = nh_trilegalisl(paths, data, tri_gaia, tri_mag, max_mag, save_file, flag_method, tri_type)

fieldnum = data.header.fieldnum;

% step 1: compute the area of a LORRI image and pull out that many stars
% from the list at random
surveyarea = data.astrom.imagew.*data.astrom.imageh.*...
    data.cal.pixsize_arcsec.^2 ./ 3600.^2; % num of pix in image * degree/pix = total deg^2

% step 2: read in the trilegal information
% determine which files are being examined
if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/ghosts/') == 1
    ghost = 1;
    old = 0;
    new = 0;
elseif strcmp(paths.datadir,'/data/symons/NH_old_data/mat/good/') == 1
    old = 1;
    ghost = 0;
    new = 0;
elseif strcmp(paths.datadir,'/data/symons/nh_data/mat/') == 1
    new = 1;
    old = 0;
    ghost = 0;
elseif strcmp(paths.datadir,'/data/symons/nh_data_lauer/mat/') == 1
    new = 0;
    old = 0;
    ghost = 0;
    lauer = 1;
end
% load the appropriate trilegal catalog files
if new == 1
    load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/trilegal/new_data/%s/isltrilegal_%02d.mat',tri_type,fieldnum));
elseif old == 1
    load(sprintf('/home/symons/isl_trilegal/%s/isltrilegal_%02d.mat',tri_type,fieldnum));
elseif lauer == 1
    load(sprintf('/home/symons/nh_ebl_pipeline/matlab/lookup/trilegal/lauer_data/%s/isltrilegal_%02d.mat',tri_type,fieldnum));
end

nfiles = numel(V);

isltot = zeros(nfiles,1);
islmasked = zeros(nfiles,1);
islghost = zeros(nfiles,1);

for jfile=1:nfiles

    % If using trilegal catalog based on UBVRI magnitudes
    if strcmp(tri_type,'ubvri')
        mag = V(jfile).mlf;
        % save trilegal mags to text file
        %     dlmwrite(['trimag',num2str(jfile),'.txt'],mag,'delimiter','\n','precision',8)

        % step 3: convert from LORRI-band mag to flux
        if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
            Fcat = data.cal.vzero .* 10.^(-V(jfile).mlf/2.5);
        elseif strcmp(flag_method, 'old') == 1
            Fcat = 3055 .* 10.^(-V(jfile).mlf/2.5);
        end

        % step 4: make a mask function
        if tri_gaia == 0
            if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
                whpl = V(jfile).V > max_mag;
            elseif strcmp(flag_method, 'old') == 1
                whpl = V(jfile).V > data.mask.maxmag;
            end
        elseif tri_gaia == 1
            whpl = V(jfile).V > tri_mag;
        end
        % If using trilegal catalog based on G magnitudes
    elseif strcmp(tri_type,'gaia')
        mag = V(jfile).G;
        % save trilegal mags to text file
        %     dlmwrite(['trimag',num2str(jfile),'.txt'],mag,'delimiter','\n','precision',8)

        % step 3: convert from LORRI-band mag to flux
        if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
            Fcat = data.cal.vzero .* 10.^(-V(jfile).G/2.5);
        elseif strcmp(flag_method, 'old') == 1
            Fcat = 3055 .* 10.^(-V(jfile).G/2.5);
        end

        % step 4: make a mask function
        if tri_gaia == 0
            if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1)
                whpl = V(jfile).G > max_mag;
            elseif strcmp(flag_method, 'old') == 1
                whpl = V(jfile).G > data.mask.maxmag;
            end
        elseif tri_gaia == 1
            whpl = V(jfile).G > tri_mag;
        end
    end

    % Select only stars with Gaia mag > 8 (used for diffuse ghost calc)
    %     diffghoststars = (V(jfile).mlf-data.cal.gaia2lorrimag) > 8;

    magcat = mag(whpl);

    % step 5: convert to surface brightness
    lIltot = 1e-26.*1e9.*data.cal.nu.*Fcat./(surveyarea .* (pi./180).^2);
    lIlcat = 1e-26.*1e9.*data.cal.nu.*Fcat(whpl)./(surveyarea .* (pi./180).^2); %this one is ISL - per radian squared
    %     lIlghost = 1e-26.*1e9.*data.cal.nu.*Fcat(diffghoststars)./(surveyarea .* (pi./180).^2); %diffuse ghost stars only

    isltot(jfile) = sum(lIltot);
    islmasked(jfile) = sum(lIlcat); %this is ISL
    %     islghost(jfile) = sum(lIlghost);%total diffuse ghost star ISL

end

% star_list = horzcat(magcat,(lIlcat.*(data.astrom.imagew.*data.astrom.imageh)));
%
% if tri_gaia == 0
%     filename = strcat(num2str(max_mag),'_tri_only.mat');
% elseif tri_gaia == 1
%     filename = strcat(num2str(tri_mag),'_tri.mat');
% end
%
% if save_file == 1
%     save(filename,'star_list');
% end

islout.isltotmean = mean(isltot);
islout.isltoterr = std(isltot);

islout.islmaskedmean = mean(islmasked); %this is ISL
islout.islmaskederr = std(islmasked);

% % Calculate total ISL from stars with Gaia mag > 8 (used for diffuse ghost)
% totislghost = mean(islghost);
% % Calculate estimate of expected summed diffuse ghost contribution based on
% % slope of ghost sb vs. star sb
% ghostsbcomp = data.ghost.diffuseslope*totislghost;
% % Compare to summed diffuse ghost sb
% ghostsbdiff = data.ghost.diffusesub - ghostsbcomp;
% % Save comparison to data
% data.ghost.diffuseisl = totislghost;
% data.ghost.diffusecomp = ghostsbcomp;
% data.ghost.diffusediff = ghostsbdiff;

% fileout = sprintf('%sisl/%s',paths.tridir,datafiles(ifile).name);
%
% save(fileout,'islout')

data.isl.tritotmean = islout.isltotmean;
data.isl.tritoterr = islout.isltoterr;
data.isl.trimean = islout.islmaskedmean; %this is used as ISL subtraction
data.isl.trierr = islout.islmaskederr;

end


