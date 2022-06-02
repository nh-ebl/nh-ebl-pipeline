function catalog_data_gaia() %gals,paths

% set variables to run manually instead of in pipeline
paths = get_paths_newest(); % CHANGE FIELD COORDS BELOW
gals = 0;

%retrieve list of gaia catalog files
cat_files = dir(fullfile(sprintf('%scatalog_files/',paths.gaiadir),'*.fit'));
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
    fprintf('Current field: %d\n', field);
    
    %get catalog file from that field
    data = fitsread(sprintf('%scatalog_files/%s',paths.gaiadir,fieldName),'asciitable');
    
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
    if not(isfolder([paths.gaiadir,'mat_files']))
        mkdir([paths.gaiadir,'mat_files'])
    end
    save(sprintf('%smat_files/field_%d_data', paths.gaiadir,field), 'RA', 'DEC', 'Gmag', 'Gmagerr');
    field_number=[field_number field];
end

%manually save field ra, dec, and number to cataloginfo mat file [MUST
%enter in CONSECUTIVE field numbers e.g. 1,2,3,4... EVEN IF field_num var
%is NOT consecutive! (that's fixed later)!}

% new fields
% field_DEC=[field_DEC -9.36002165 4.792341167 -22.82961047 25.90105487 4.828636396 -22.80892287 25.9444605 12.26639126];
% field_RA=[field_RA 258.746143 220.7927148 191.3474659 259.7907857 220.7621693 191.3286772 259.8207519 235.187082];

% lauer fields
% field_DEC=[field_DEC -17.7780 -41.6360 -50.7529 -0.5183 0.2915 36.2331 35.2979];
% field_RA=[field_RA 1.7790 348.0611 33.4069 358.2428 0.8066 224.9875 226.4865];

% newest fields
field_DEC=[field_DEC -0.917613 38.383623 -14.025399 -13.8919845 -19.8081555 -23.432433 11.0600625 59.8891035 34.180979 33.354217 32.6292905 21.3322165 12.3094285 33.7912115 -53.863063 14.898267 -48.918064 -40.4043165 0.6843575 -42.9560245 -44.888367 -34.836905 -59.4490655];
field_RA=[field_RA 22.536525 253.81149 344.507784 345.007431 315.679527 342.440189 209.0270525 234.477303 228.930197 229.9648425 230.807135 242.346948 249.9032945 351.786695 331.469115 230.3917645 312.2470555 295.680608 247.6972585 352.764887 49.6211205 296.866018 333.41013];

%complicated reorder to deal with skipped field numbers and linux file read order weirdness
field_DEC_reorder = zeros(size(field_number));
field_RA_reorder = zeros(size(field_number));
field_number_sorted = sort(field_number);
for i = 1:length(field_DEC)
    field_DEC_reorder(i) = field_DEC( field_number(i) == field_number_sorted );
    field_RA_reorder(i) = field_RA( field_number(i) == field_number_sorted );
end
field_DEC = field_DEC_reorder;
field_RA = field_RA_reorder;

save(sprintf('%s/mat_files/cataloginfo',paths.gaiadir), 'field_number','field_RA','field_DEC');
end

