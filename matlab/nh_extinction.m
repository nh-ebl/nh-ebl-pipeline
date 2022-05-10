function data = nh_extinction(data, paths)

%check if ext maps are already made, if not call python script that makes
%them (isfile returns 1 if file exists)
if( ~(isfile(sprintf('%sext_%s.fits',paths.extdir,data.header.timestamp))))
    %if at least one of the files isn't there, call python script to make
    %them all
    disp('No ext file found, retrieving new ext file.')
    pydir = '/home/symons/nh_ebl_pipeline/py/extinction/'; %where the python script is
    pyfile = 'get_extinction.py'; %name of the python file to run
    imagefile = [paths.imagedir,'regist_',data.header.rawfile]; %fits file to send to the python script
    %write the imagefile path in a text file to where the python script is
    fileID = fopen([pydir,'imagefile.txt'],'w'); %open the file to be written
    fprintf(fileID,'%s',imagefile); %write imagefile path
    fclose(fileID); %close file
    %call python script
    system(['python ',pydir,pyfile]);
    if not(isfolder([paths.extdir]))
        mkdir([paths.extdir])
    end
    %get the ext files and move them to their data directory
    movefile([pydir,'ext_',data.header.timestamp,'.fits'], [paths.extdir,'ext_',data.header.timestamp,'.fits']); %move from pydir to paths.extdir
end

% Load in extinction map
extmap = fitsread(sprintf('%sext_%s.fits',...
    paths.extdir,data.header.timestamp));

% Calculate flux or intensity ratio
flux_rat = 10.^(-0.4.*extmap);
mean_flux_rat = mean(mean(flux_rat));

% Save to data
data.ext.extmap = extmap;
data.ext.flux_rat = flux_rat;
data.ext.mean_flux_rat = mean_flux_rat;
end

