function data=nh_calcisl(data, paths, params, use_gaia, tri_gaia, tri_mag, wing_mag, max_mag, save_file, flag_method,errflag_psf,tri_type)

% load('run_params.mat','params')

% Calculate ISL from PSF wings for stars < X mag (using Gaia or USNOB1)
usnoisl = nh_usnoisl(data, paths, use_gaia, wing_mag, save_file, flag_method, errflag_psf);

if params.err_psf == 1
    data.(params.err_str).isl.usnotot = usnoisl.isltot;
    data.(params.err_str).isl.usnototim = usnoisl.totimage;
    data.(params.err_str).isl.usnowing = usnoisl.islwing; %this one gets used
    data.(params.err_str).isl.usnowingim = usnoisl.wingimage;
    data.(params.err_str).isl.usnofaint = usnoisl.islfaint;
else
    data.isl.usnotot = usnoisl.isltot;
    data.isl.usnototim = usnoisl.totimage;
    data.isl.usnowing = usnoisl.islwing; %this one gets used
    data.isl.usnowingim = usnoisl.wingimage;
    data.isl.usnofaint = usnoisl.islfaint;
end

%   gsciiisl = nh_gsciiisl(data);
%
%   data.isl.gsciitot = gsciiisl.isltot;
%   data.isl.gsciitotim = gsciiisl.totimage;
%   data.isl.gsciiwing = gsciiisl.islwing;
%   data.isl.gsciiwingim = gsciiisl.wingimage;
%   data.isl.gsciifaint = gsciiisl.islfaint;

if tri_gaia == 1
    
    % Calculate ISL from gaia for mask mag to tri mag and ISL from trilegal
    % for mag > tri mag
    data = nh_gaiaisl(data, paths, tri_mag, max_mag, save_file);
    data = nh_trilegalisl(paths, data, tri_gaia, tri_mag, max_mag, save_file);
    
elseif tri_gaia == 0
    
    % Calculate ISL from trilegal for mag > mask mag
    data = nh_trilegalisl(paths, data, tri_gaia, tri_mag, max_mag, save_file, flag_method,tri_type);
    
    %     triout=nh_add_trilegalisl(data, paths);
    %     data.isl.tritotmean = triout.isltotmean;
    %     data.isl.tritoterr = triout.isltoterr;
    %     data.isl.trimean = triout.islmaskedmean;
    %     data.isl.trierr = triout.islmaskederr;
end

end