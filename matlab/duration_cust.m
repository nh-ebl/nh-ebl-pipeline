% MATLAB script to extract exposure time of NH images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments added by C.Nguyen
% Revision 11/04/15: comments added, file paths updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
% data = path to data directory
% ndate = number of files in data directory
% date_array = array of the dates of when the data was taken
% file_array = array of the file names
% duration = array of exposure time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function duration(phase,path)

% open the data directory 
data = dir(sprintf('%s%s/data',path,phase)); 
ndata=numel(data);

% define arrays
date_array=[]; 
file_array=[];
durations=[];
k=1;

% 1st loop over every .fit in the data directory
for idate= 1:ndata;
    % pass the .fit name to filename
    filename = sprintf('%s%s/data/%s/*fit',path,phase,data(idate).name);
    files = dir(filename)     ;
    
    % 2nd loop over the header keywords
    % to find the exposure time
    for ifile = 1:numel(files)
       date_array=[date_array data(idate)];
       file_array=[file_array files(ifile)] ;
       
      
       info= fitsinfo(sprintf('%s%s/data/%s/%s',path,phase,data(idate).name,files(ifile).name));
       exptime= info.PrimaryData.Keywords(22,2);
       durations = [durations exptime];
       k=k+1;
    end 
   
end 

duration_array= cell2mat(durations);

% save the extracted information as a matrix for MATLAB to deal with later
save(sprintf('%smatrix/NH_%s_time.mat',path,phase), 'date_array', 'file_array', 'duration_array');

return 
 
     
        
    
    