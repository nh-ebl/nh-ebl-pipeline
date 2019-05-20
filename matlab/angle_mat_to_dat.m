%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C.Nguyen
% MATLAB script to export MATLAB matrix into .dat file if needed
 function angle_mat_to_dat(phase)

launch_date = 2453755.291667;
au = 1.496*10^8;

filename = fopen(sprintf('../../data/NH/txt/%s_sun_NH.txt',phase), 'w');
fprintf(filename, '%%t t_from_launch x y z d\n');

load (sprintf('../../data/NH/matrix/%s_d_NH_sun.mat',phase));
date = date_array;
file = file_array;
xd = xd_sun_NH_array;
yd = yd_sun_NH_array;
zd = zd_sun_NH_array;

[Nrows,Ncols] = size(xd);

    for i=1:Ncols
        
        date_from_launch = date(i) - launch_date ;
        
        d(i) = sqrt(xd(i)^2 + yd(i)^2 + zd(i)^2)/au;
        xd(i) = xd(i)/au;
        yd(i) = yd(i)/au;
        zd(i) = zd(i)/au;
        
        fprintf(filename,'%d %d %d %d %d %d \n', date(i), date_from_launch, xd(i), yd(i), zd(i), d(i));
    end

fclose(filename);

filename = fopen(sprintf('../../data/NH/txt/%s_angle_NH_Earth_sun.txt',phase), 'w');
fprintf(filename, '%%t t_from_launch toEarth toSun\n');

load (sprintf('../../data/NH/matrix/%s_angle_NH_Earth_sun.mat',phase));

date = date_array;
file = file_array;
Earth = angle_NH_Earth_array;
Sun = angle_NH_sun_array;

[Nrows,Ncols] = size(xd);

    for i=1:Ncols
        
        date_from_launch = date(i) - launch_date ;
        
        
        fprintf(filename,'%d %d %d %d\n', date(i), date_from_launch, Earth(i), Sun(i));
    end


fclose(filename);
return

