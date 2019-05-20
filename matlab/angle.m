%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C.Nguyen
% MATLAB Script to extract information about target and spacecraft
% distance as a function of time for NH LORRI data. 
% The data sets are selected by launch_selection.m and launch_selection.m
% 11/06/2015: ver 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
% data = path to data directory
% ndate = number of files in data directory
% date_array = array of the dates of when the data was taken JD
%
% file_array = array of the file names
%
% unit of distances [km]
% xd_NH_target_array = array of the x-distance from the spacecraft to target
% yd_NH_target_array = array of the y-distance from the spacecraft to target
% zd_NH_target_array = array of the z-distance from the spacecraft to target
%
% xd_sun_target_array = array of the x-distance from the sun to target
% yd_sun_target_array = array of the y-distance from the sun to target
% zd_sun_target_array = array of the z-distance from the sun to target
%
% xd_NH_sun_array = array of the x-distance from the sun to the spacecraft
% yd_NH_sun_array = array of the x-distance from the sun to the spacecraft
% zd_NH_sun_array = array of the x-distance from the sun to the spacecraft
%
% target_name_array = array of target names. Names are in SPICE catalog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%import matlab.io.*;

function angle(phase)
import matlab.io.*;

% open the data directory 
data = dir(sprintf('../../data/NH/%s/selected_data', phase)); 
ndata=numel(data);

% define arrays
date_array=[]; 
file_array=[];

angle_NH_Earth_array=[];
angle_NH_sun_array=[];

xd_sun_NH_array=[];
yd_sun_NH_array=[];
zd_sun_NH_array=[];

k=1;

% 1st loop over every .fit in the data directory
    % pass the .fit name to filename
    filename = sprintf('../../data/NH/%s/selected_data/*fit',phase) ;    
    files = dir(filename);
    
    % 2nd loop over the header keywords
    % to find the distance
    for ifile = 1:numel(files);
       file_array=[file_array files(ifile)] ;
       
       % open fits file
       file_to_read= fits.openFile(sprintf('../../data/NH/%s/selected_data/%s',phase, files(ifile).name));
       
       % search for values using keywords
       date= fits.readKey(file_to_read, 'SPCUTCJD');
       date_split = strsplit(date,' ');
       date_no_jd = str2double(date_split(2));
       date_array = [date_array date_no_jd];
       
       % extract angle from NH to Earth and sun
       angle_NH_Earth = fits.readKey(file_to_read, 'SOL_ELON');
       angle_NH_Earth = str2num(angle_NH_Earth);
       angle_NH_Earth_array = [angle_NH_Earth_array angle_NH_Earth];
       angle_NH_sun = fits.readKey(file_to_read, 'EAR_ELON');
       angle_NH_sun = str2num(angle_NH_sun);
       angle_NH_sun_array = [angle_NH_sun_array angle_NH_sun];
             
       
       % extract distance from NH to sun
       x_sun_NH = fits.readKey(file_to_read, 'SPCSSCX');
       x_sun_NH = str2num(x_sun_NH);
       xd_sun_NH_array = [xd_sun_NH_array x_sun_NH];
       y_sun_NH = fits.readKey(file_to_read, 'SPCSSCY');
       y_sun_NH = str2num(y_sun_NH);
       yd_sun_NH_array = [yd_sun_NH_array y_sun_NH];
       z_sun_NH = fits.readKey(file_to_read, 'SPCSSCZ');
       z_sun_NH = str2num(z_sun_NH);
       zd_sun_NH_array = [zd_sun_NH_array z_sun_NH];
       
       fits.closeFile(file_to_read);
       
       k=k+1;
    end 

% save the extracted information as matrices for MATLAB to deal with later
save(sprintf('../../data/NH/matrix/%s_d_NH_sun.mat',phase), 'date_array', 'file_array', 'xd_sun_NH_array', 'yd_sun_NH_array', 'zd_sun_NH_array');
save(sprintf('../../data/NH/matrix/%s_angle_NH_Earth_sun.mat',phase), 'date_array', 'file_array', 'angle_NH_Earth_array', 'angle_NH_sun_array');

return