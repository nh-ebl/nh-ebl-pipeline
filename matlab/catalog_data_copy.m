function catalog_data_copy()

n = 11; %number of fields to look at (stored as folders name 1 to n

%column numbers for ra, dec, and mag values in USNOB-1 ascii files
% there s/g values are + 3 colunmns from each mag val.
ra = 2;
dec = 3;
b1 = 4;
b2 = 10;
r1 = 7;
r2 = 13;
i2 = 16;
sg = 2;

field_DEC = [];
field_RA = [];
field_number = [];
%loop through each field
for i=1:n
	field = i;
    fprintf('Current field: %i\n', field); 
	
	%get all catalog files from that field
	data = dir(sprintf('C:\\nh_data\\usnob1\\catalog_files\\%d\\*field*',i));
	ndata = numel(data);
    fprintf('Number of subfields: %i\n', ndata); 
	
	%initialize array values
	RA = [];
	DEC = [];
	B1mag = [];
	B2mag = [];
	R1mag = [];
	R2mag = [];
	I2mag = [];
	B1sg = [];
	B2sg = [];
	R1sg = [];
	R2sg = [];
	I2sg = [];
	
	%loop through each feild for current field
	for j=1:ndata
        fprintf('Current subfield: %i\n', j); 
		%open each file
		filename = data(j).name;
		fileID = fopen(sprintf('C:\\nh_data\\usnob1\\catalog_files\\%d\\%s', i, filename));
		
		%read first line
		tline = fgets(fileID);
	
		%we want to ignore the lines until the table starts 
		collect = 0;
		
		%loop through all lines in current file
		while tline ~= -1
			C = strsplit(tline);
			%C = cell2mat(C)
			%if we have reached table, retrieve relevent data
			%otherwise check if we have reach the line before the 
			%table starts. This line starts with '#3'
            %1 + the column number since there is a whitespce charecter 
            %when the line is split
			if collect == 1
				RA = [RA str2double(cell2mat(C(ra + 1)))];
				DEC = [DEC str2double(cell2mat(C(dec + 1)))];
				B1mag = [B1mag str2double(cell2mat(C(b1 + 1)))];
				B2mag = [B2mag str2double(cell2mat(C(b2 + 1)))];
				R1mag = [R1mag str2double(cell2mat(C(r1 + 1)))];
				R2mag = [R2mag str2double(cell2mat(C(r2 + 1)))];
				I2mag = [I2mag str2double(cell2mat(C(i2 + 1)))];
				B1sg = [B1sg str2double(cell2mat(C(b1 + sg + 1)))];
				B2sg = [B2sg str2double(cell2mat(C(b2 + sg + 1)))];
				R1sg = [R1sg str2double(cell2mat(C(r1 + sg + 1)))];
				R2sg = [R2sg str2double(cell2mat(C(r2 + sg + 1)))];
				I2sg = [I2sg str2double(cell2mat(C(i2 + sg + 1)))];
			else 
				if strcmp(C(1), '|null') == 1
					collect = 1; %set collect to 1
				end
			end
			tline = fgets(fileID);
		end
		fclose(fileID);		
    end	
	%save data from field as matlab matrices
	save(sprintf('C:\\nh_data\\usnob1\\mat_files\\field_%d_data', i), 'RA', 'DEC', ...
        'B1mag', 'B2mag', 'R1mag', 'R2mag', 'I2mag', ...
        'B1sg', 'B2sg', 'R1sg', 'R2sg', 'I2sg');
	%[RA' DEC' B1mag' B2mag' R1mag' R2mag'];
    field_number=[field_number i];
end
field_DEC=[field_DEC -9.36002165 4.792341167 -22.82961047 25.90105487 4.828636396 -22.80892287 25.9444605 12.26639126 6.735450267 43.2828838 6.793556222];
field_RA=[field_RA 258.746143 220.7927148 191.3474659 259.7907857 220.7621693 191.3286772 259.8207519 235.187082 136.1413223 157.6944323 136.1518358];
save(sprintf('C:\\nh_data\\usnob1\\mat_files\\cataloginfo'), 'field_number','field_RA','field_DEC');
	
