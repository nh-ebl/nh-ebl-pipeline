function texp_vs_t_ext(phase,path)

import matlab.io.*;

data = dir(sprintf('%s%s/data',path,phase)); 
ndata=numel(data);

t_exp = [];
date_array =[];

k=1;

for idate= 1:ndata;
    % pass the .fit name to filename
    filename = sprintf('%s%s/data/%s/*fit',path,phase,data(idate).name);
    files = dir(filename)     ;
    
    % 2nd loop over the header keywords
    % to find the exposure time
    for ifile = 1:numel(files);
        
       file_to_read= fits.openFile(sprintf('%s%s/data/%s/%s',path, phase,data(idate).name, files(ifile).name));
      
       date= fits.readKey(file_to_read, 'SPCUTCJD');
       date_split = strsplit(date,' ');
       date_no_jd = str2double(date_split(2));
       date_array = [date_array date_no_jd];
       
       time= fits.readKey(file_to_read, 'EXPTIME');
       time_num = str2double(time);
       t_exp = [t_exp time_num];
       
       fits.closeFile(file_to_read);
       
       k = k+1;
       
    end 
   
end 

save(sprintf('%smatrix/%s_texp_t.mat',path,phase), 'date_array', 't_exp');

return