function data = nh_scattering(data, paths)

%check if planck maps are already made, if not call python script that makes
%them (isfile returns 1 if file exists)
if( ~(isfile(sprintf('%s%s.mat',paths.scatteringdir,data.header.timestamp))))
    %if at least one of the files isn't there, call python script to make
    %them all
    disp('No scattering file found, creating new scattering file.')
    pydir = '/home/symons/nh_ebl_pipeline/py/scattering/'; %where the python script is
    pyfile = 'get_scattering_image.py'; %name of the python file to run
    imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
    %write the imagefile path in a text file to where the python script is
    fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
    fprintf(fileID,'%s',imagefile); %write imagefile path
    fclose(fileID); %close file
    %call python script
    system(['python ',pydir,pyfile]);
    %get the scattering files and move them to their data directory
    movefile([pydir,data.header.timestamp,'.mat'], [paths.scatteringdir,data.header.timestamp,'.mat']); %move from pydir to paths.scatteringdir
end

lgmap_nano = load([paths.scatteringdir,data.header.timestamp,'.mat']);

end