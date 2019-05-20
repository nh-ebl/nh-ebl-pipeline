function mean_flux_mat_to_dat(phase)
load (sprintf('../../data/NH/matrix/%s_mean_flux.mat',phase));

[Nrows,Ncols] = size(mean_array);

filename = fopen(sprintf('../../data/NH/txt/mean_flux_non_physical.txt',phase), 'w');
fprintf(filename, '%%t t_from_launch d mean stdev\n');

for i=1:Ncols
    fprintf(filename,'%d %d %d %d %d \n', date_array(i), date_from_launch(i), d(i), mean_array(i), sigma_array(i));
end

fclose(filename);

return

