%test file for retriving RA/Dec from .fit data

funtion test_ra_da(phase)
	import matlab.io.*;

	RA = [];
	DEC = [];
	file = [];
	date = [];

	filename = sprintf('../../data/NH/%s/selected_data/*fit', phase);
	files = dir(filename); %??

	for i = 1:numel(files);
	
		file = [file files(i)];

		file_to_read = fits.openFile(sprintf('../../data/NH/%s/selected_data/%s,phase,files(i).name));

		currdate = fits.readKey(file_to_read, 'SPCUTJD');
		date_split = strsplit(currdata,' ');
		date_jd = str2double(date_split(2));
		date = [date date_jd];
	
		currRA = fits.readKey(file_to_read, 'CRVAL1');
		RA_num = str2double(currRA);
		RA = [RA currRA];
	
		currDEC = fits.readKet(file_to_read, 'CRVAL2');
		DEC_num = string2double(currDEC);
		DEC = [DEC DEC_num];

		fits.closeFile(file_to_read);

	end
	
return
	

 

