function catalog_data_gaia_mc() %gals,paths

% set variables to run manually instead of in pipeline
% paths = get_paths_lauer();
% paths = get_paths_new();
% paths = get_paths_old();
paths = get_paths_newest();

gals = 0;
mc_num = 100;

% Set the random seed
s = RandStream('mlfg6331_64','Seed',1776);

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
    for m = 1:mc_num
        fprintf('Current mc num: %i\n',m);
        fieldName = cat_filesTemp{1,i}; %get field fit name
        field = str2num(cat_filesTemp{1,i}(6:strfind(cat_filesTemp{1,i},'.fit')-1));  %get field fit number from fit name
        fprintf('Current field: %i\n', field);

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

            % Figure out how many entries to be randomly chosen
            tot_len = length(name_gal); % Total number of entries
            test_num = round(tot_len*0.2771); % We want 27.71% (purity) of total - see Bailer-Jones et al.
            get_ind = randsample(s,tot_len,test_num) % Randomly generate selected indices from seed

            %Remove randomly chosen galaxies from list of all objects
            for j=1:length(get_ind)
                k = contains(name_direct,name_gal(get_ind(j))); % Find the name_direct that contains a name_gal
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
        if not(isfolder([paths.gaiadir,'mat_files/field_',num2str(field),'_mc']))
            mkdir([paths.gaiadir,'mat_files/field_',num2str(field),'_mc'])
        end
        save(sprintf('%smat_files/field_%d_mc/%i', paths.gaiadir,field,m), 'RA', 'DEC', 'Gmag', 'Gmagerr');
        %     field_number=[field_number field];
    end

end

%manually save field ra, dec, and number to cataloginfo mat file

% new fields
% field_DEC=[field_DEC -9.36002165 4.792341167 -22.82961047 25.90105487 4.828636396 -22.80892287 25.9444605 12.26639126];
% field_RA=[field_RA 258.746143 220.7927148 191.3474659 259.7907857 220.7621693 191.3286772 259.8207519 235.187082];

% lauer fields
% field_DEC=[field_DEC -17.7780 -41.6360 -50.7529 -0.5183 0.2915 36.2331 35.2979];
% field_RA=[field_RA 1.7790 348.0611 33.4069 358.2428 0.8066 224.9875 226.4865];
%
% save(sprintf('%s/mat_files/cataloginfo',paths.gaiadir), 'field_number','field_RA','field_DEC');
end

