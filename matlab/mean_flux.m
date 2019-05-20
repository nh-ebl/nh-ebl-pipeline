function mean_flux(phase)
distance_for_plot('regist');
load (sprintf('../../data/NH/matrix/%s_distance_date.mat',phase));

filename = sprintf('../../data/NH/%s/masked/*fit',phase) ;  

files = dir(filename);

mean_array = [];
stdev_array=[];
sigma_array=[];
N_pixel=[];

lambda_pivot = 6076; % [A]
dlambda = 5000;
convert = 10.0^(6); % convert Wcm-2sr-1 to nWm-2sr-1
rsun = 266400;

for ifile = 1:numel(files);
    file_to_read= sprintf('../../data/NH/%s/masked/%s',phase,files(ifile).name);
    data_raw = single(fitsread(file_to_read));
    data = reshape(data_raw, [numel(data_raw),1]);
    
    % find mean pixel counts
    mean_array(ifile) = nanmean(data)*convert*dlambda/rsun;
    
    stdev_array(ifile) = nanstd(data)*convert*dlambda/rsun;
    
    N_pixel(ifile) = sum(~isnan(data));
    sigma_array(ifile) = sqrt((stdev_array(ifile))^2.0/sqrt(N_pixel(ifile)));
    
end

save(sprintf('../../data/NH/matrix/%s_mean_flux.mat',phase),'date_array','date_from_launch', 'd', 'mean_array','sigma_array')
return

function distance_for_plot(phase)
distance(phase);

load (sprintf('../../data/NH/matrix/%s_sun_NH.mat',phase));

launch_date = 2453755.291667;
au = 1.496*10^8;


xd = xd_sun_NH_array;
yd = yd_sun_NH_array;
zd = zd_sun_NH_array;
date = date_array;
date_from_launch=[];
d=[];

[Nrows,Ncols] = size(xd); 

for i=1:Ncols
   
        date_from_launch(i) = date(i) - launch_date ;
        
        d(i) = sqrt(xd(i)^2 + yd(i)^2 + zd(i)^2)/au;
        xd(i) = xd(i)/au;
        yd(i) = yd(i)/au;
        zd(i) = zd(i)/au;
end
save(sprintf('../../data/NH/matrix/%s_distance_date.mat',phase), 'date_array','date_from_launch','d')
return
