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

function distance(phase,path)
import matlab.io.*;

% open the data directory 
data = dir(sprintf('%s%s/selected_data',path,phase)); 
ndata=numel(data);

% define arrays
date_array=[]; 
file_array=[];

xd_NH_target_array=[];
yd_NH_target_array=[];
zd_NH_target_array=[];

xd_sun_target_array=[];
yd_sun_target_array=[];
zd_sun_target_array=[];

xd_sun_NH_array=[];
yd_sun_NH_array=[];
zd_sun_NH_array=[];

target_name_array=[];

k=1;

% 1st loop over every .fit in the data directory
    % pass the .fit name to filename
    filename = sprintf('%s%s/selected_data/*fit',path,phase) ;    
    files = dir(filename);
    
    % 2nd loop over the header keywords
    % to find the distance
    for ifile = 1:numel(files);
       file_array=[file_array files(ifile)] ;
       
       % open fits file
       file_to_read= fits.openFile(sprintf('%s%s/selected_data/%s',path,phase,files(ifile).name));
       
       % search for values using keywords
       date= fits.readKey(file_to_read, 'SPCUTCJD');
       date_split = strsplit(date,' ');
       date_no_jd = str2double(date_split(2));
       date_array = [date_array date_no_jd];
       
       % extract distance from target to NH
       x_NH_target = fits.readKey(file_to_read, 'SPCTSCX');
       x_NH_target = str2num(x_NH_target);
       xd_NH_target_array = [xd_NH_target_array x_NH_target];
       y_NH_target = fits.readKey(file_to_read, 'SPCTSCY');
       y_NH_target = str2num(y_NH_target);
       yd_NH_target_array = [yd_NH_target_array y_NH_target];
       z_NH_target = fits.readKey(file_to_read, 'SPCTSCZ');
       z_NH_target = str2num(z_NH_target);
       zd_NH_target_array = [zd_NH_target_array z_NH_target];
             
       % extract distance from target to sun
       x_sun_target = fits.readKey(file_to_read, 'SPCTSOX');
       x_sun_target = str2num(x_sun_target);
       xd_sun_target_array = [xd_sun_target_array x_sun_target];
       y_sun_target = fits.readKey(file_to_read, 'SPCTSOY');
       y_sun_target = str2num(y_sun_target);
       yd_sun_target_array = [yd_sun_target_array y_sun_target];
       z_sun_target = fits.readKey(file_to_read, 'SPCTSOZ');
       z_sun_target = str2num(z_sun_target);
       zd_sun_target_array = [zd_sun_target_array z_sun_target];
       
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
       
       % extract target names
       target_name = fits.readKey(file_to_read, 'SPCTCB');
       target_name_array = strvcat(target_name_array, target_name);
       
       fits.closeFile(file_to_read);
       
       k=k+1;
    end 

% save the extracted information as matrices for MATLAB to deal with later
save(sprintf('%smatrix/%s_NH_target.mat',path,phase), 'date_array', 'file_array', 'xd_NH_target_array', 'yd_NH_target_array', 'zd_NH_target_array');
save(sprintf('%smatrix/%s_sun_target.mat',path,phase), 'date_array', 'file_array', 'xd_sun_target_array', 'yd_sun_target_array', 'zd_sun_target_array');
save(sprintf('%smatrix/%s_sun_NH.mat',path,phase), 'date_array', 'file_array', 'xd_sun_NH_array', 'yd_sun_NH_array', 'zd_sun_NH_array');
save(sprintf('%smatrix/%s_target_list.mat',path,phase), 'target_name_array');
return