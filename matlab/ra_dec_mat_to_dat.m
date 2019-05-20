function ra_dec_mat_to_dat(phase,path)

load (sprintf('%smatrix/%s_ra_dec.mat',path,phase));

filename = fopen(sprintf('%stxt/%s_RA_DEC.txt',path,phase), 'w');
fprintf(filename, '%%t t_from_launch RA DEC\n');


date = date_array;
lat=[];
long=[];

e = deg2rad(23.4);
launch_date = 2453755.291667;

RA = 0.0;
DEC = 0.0;
date = 0.0;

[Nrows,Ncols] = size(date_array);

for i=1:Ncols
    date = date_array(i) - launch_date;
    fprintf(filename, '%s %s %s %s\n', date_array(i), date, RA_array(i), DEC_array(i));
end

fclose(filename);

return
