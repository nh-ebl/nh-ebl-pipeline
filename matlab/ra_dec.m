%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C.Nguyen
% MATLAB Script to extract information about target RA & DEC
% the output matrix will then be fed into an astrometric registration code
% 02/06/2017: ver 2.0
% update: added JD, T_exposure, and target field
% update: header keywords for the target field RA/DEC changed from that of
% the reference pixels to the instrument pointing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
% data = path to data directory
% ndate = number of files in data directory
% RA_array = array of the image center RA
% DEC_array = array of the image center DEC
% file_array = array of filename, needed for astrometric registration later
% exp_array = array of exposure time
% date_array = array of JD
% target_array = array of targets of observations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%import matlab.io.*;

function ra_dec(phase,path)
import matlab.io.*;

% open the data directory 
data = dir(sprintf('%s%s/selected_data',path,phase)) ;
ndata=numel(data);

% define arrays
RA_array=[]; 
DEC_array=[];
file_array=[];
file_name_array=[];
date_array=[];
exp_array=[]; % exposure time
target_array=[]; % target name

k=1;

% 1st loop over every .fit in the data directory
    % pass the .fit name to filename
    filename = sprintf('%s%s/selected_data/*fit',path,phase); 
    files = dir(filename);
    
    for ifile = 1:numel(files);
        
       file_array = [file_array files(ifile)];
       file_name = files(ifile).name;
       file_name_array{end+1} = file_name;
       
       % open fits file
       file_to_read= fits.openFile(sprintf('%s%s/selected_data/%s',path,phase,files(ifile).name));
       
       % search for values using keywords
       date= fits.readKey(file_to_read, 'SPCUTCJD');
       date_split = strsplit(date,' ');
       date_no_jd = str2double(date_split(2));
       date_array = [date_array date_no_jd];
       
       % search for values using keywords
       RA= fits.readKey(file_to_read, 'SPCBRRA');
       RA_number = str2double(RA);
       RA_array = [RA_array RA_number];
       
       DEC = fits.readKey(file_to_read, 'SPCBRDEC')  ;
       DEC_number = str2double(DEC);
       DEC_array = [DEC_array DEC_number];
       
       exp = fits.readKey(file_to_read, 'EXPTIME')  ; 
       exp_number = str2double(exp);
       exp_array = [exp_array exp_number];
 
       target = fits.readKey(file_to_read, 'SPCTCB');  
       target_array{end+1} = target;
       
       fits.closeFile(file_to_read);
       
       k=k+1;
    end 

% save the extracted information as matrices for MATLAB to deal with later
save(sprintf('%smatrix/%s_ra_dec.mat',path,phase), ...
    'file_name_array', 'date_array', 'RA_array', 'DEC_array', 'exp_array',...
    'target_array'); %'target_RA_array', 'target_DEC_array');

% export data matrix into .txt file
filename = fopen(sprintf('%stxt/%s_RA_DEC_JD.txt',path,phase), 'w');
load (sprintf('%smatrix/%s_ra_dec.mat',path,phase));

fprintf(filename,'%date\tt_exp\tRA\tDEC\n');
i = 1;
date = date_array;
file = file_name_array;

[Nrows,Ncols] = size(RA_array);

    for i=1:Ncols
        fprintf(filename,'%s\t%f\t%f\t%f\t%f\t%s\n', ...
            file{i}, date(i), exp_array(i), RA_array(i), DEC_array(i), ...
            target_array{i});
    end
return