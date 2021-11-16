function catalog_data_gaia_wide(gals,paths)

%retrieve list of gaia catalog files
cat_files = dir(fullfile(sprintf('%swide_files/',paths.gaiadir),'*.fit'));
cat_filesTemp = struct2cell(cat_files); %can't get fields to iterate well, cells do at least

%determine number of catalog files (fields to look at)
n = numel(cat_files);

field_DEC = [];
field_RA = [];
field_number = [];

%loop through each field
for i=1:n
    fieldName = cat_filesTemp{1,i}; %get field fit name
    field = str2num(cat_filesTemp{1,i}(6:strfind(cat_filesTemp{1,i},'.fit')-1));  %get field fit number from fit name
    fprintf('Current field: %i\n', field);
    
    %get catalog file from that field
    data = fitsread(sprintf('%swide_files/%s',paths.gaiadir,fieldName),'asciitable');
    
    %save catalog data
    name = data{1,4}; %DR2Name is 4th
    RA = data{1,5}; %RA_ICRS is 5th
    DEC = data{1,7};  %DE_ICRS is 7th
    raerr = data{1,6}; %e_RA_ISCRS is 6th
    decerr = data{1,8}; %e_DE_ICRS is 8th
    Gmag = data{1,9}; %gmag is 9th
    Gmagerr = data{1,10}; %gmag err is 10th
    
    %remove galaxies if desired
    if gals == 0
        %get galaxy catalog file from that field
        tableOpts = detectImportOptions(sprintf('%sgal_files/gaiagal_field_%d.txt',paths.gaiadir,field),'Delimiter','\t','HeaderLines',46); % Detect options for the table, specify delimiter and header lines
        tableOpts = setvartype(tableOpts, 'GDR2', 'char'); % Makes galaxy ID #'s character arrays
        warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames'); % Useless warning off
        gal_data = readtable(sprintf('%sgal_files/gaiagal_field_%d.txt',paths.gaiadir,field), tableOpts); % Read the file as a table
        name_gal = table2cell(gal_data(3:size(gal_data,1),1)); % Get the names from the table
        
        %strip GAIA DR2 from each object name
        k = cell2mat(strfind(name,'DR2 '));
        name_direct = cell(length(k),1); %prep b/c cell
        for j=1:length(k)
            name_direct{j} =  name{j}(k(j)+4:end); %gets rid of GAIA DR2 bit at front
        end
        
        %Remove galaxies from list of all objects
        for j=1:length(name_gal)
            k = contains(name_direct,name_gal(j)); % Find the name_direct that contains a name_gal
            name_direct = name_direct(~k); % Remove the name_direct that contains a name_gal
            name = name(~k); %remove
            RA = RA(~k); %remove
            DEC = DEC(~k); %remove
            raerr = raerr(~k); %remove
            decerr = decerr(~k); %remove
            Gmag = Gmag(~k); %remove
            Gmagerr = Gmagerr(~k); %remove
        end
        
    end
    
    %save data from field as matlab matrices
    save(sprintf('%smat_files/field_%d_data_wide', paths.gaiadir,field), 'RA', 'DEC', 'Gmag', 'Gmagerr');
    field_number=[field_number field];
end

%manually save field ra, dec, and number to cataloginfo mat file
% field_DEC=[field_DEC -9.36002165 4.792341167 -22.82961047 25.90105487 4.828636396 -22.80892287 25.9444605 12.26639126];
% field_RA=[field_RA 258.746143 220.7927148 191.3474659 259.7907857 220.7621693 191.3286772 259.8207519 235.187082];
% save(sprintf('%s/mat_files/cataloginfo',paths.gaiadir), 'field_number','field_RA','field_DEC');
end

