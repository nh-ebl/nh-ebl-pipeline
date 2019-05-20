function JD_date(phase,path)
import matlab.io.*;


% pass the .fit name to filename
 filename = sprintf('%s%s/selected_data/*fit',path,phase) ;  
%filename = sprintf('%s/*fit', phase);
files = dir(filename);
date_array=[];
file_array=[];
% read JD values
 for ifile = 1:numel(files);
 % for ifile = 1:227;
       file_array=[file_array files(ifile)] ;
       
       % open fits file
       file_to_read= fits.openFile(sprintf('%s%s/selected_data/%s',path,phase,files(ifile).name));
       %file_to_read= fits.openFile(sprintf('%s/%s',phase,files(ifile).name));
       
       % search for values using keywords
       date= fits.readKey(file_to_read, 'SPCUTCJD');
       date_split = strsplit(date,' ');
       date_no_jd = str2double(date_split(2));
       date_array = [date_array date_no_jd];
       
       fits.closeFile(file_to_read);
       
 end 

 save(sprintf('%smatrix/%s_selected_JD_name.mat',path,phase), 'date_array', 'file_array')
%save(sprintf('%s_selected_JD_name.mat',phase), 'date_array', 'file_array')

end

