clear all
close all

%get paths for new data files or old data files
paths = get_paths_lauer();
datadir = paths.datadir;

%load directory of saved mat files for all images
datafiles = dir(sprintf('%s*.mat',datadir));

% Preallocate
exptime = zeros(size(datafiles));
mean_dn_s = zeros(size(datafiles));
mean_cal = zeros(size(datafiles));
fieldnum = zeros(size(datafiles));
mean_dn = zeros(size(datafiles));
im_num = zeros(size(datafiles));
header_bias = zeros(size(datafiles));
calc_bias = zeros(size(datafiles));
bias_diff = zeros(size(datafiles));
bias_type = zeros(size(datafiles));

fields = [1,2,3,4,5,6,7];

field_mean_dn_s = zeros(size(fields));
field_mean_dn = zeros(size(fields));
field_mean_cal = zeros(size(fields));
field_dn_s = zeros(size(fields));
field_dn = zeros(size(fields));
field_cal = zeros(size(fields));


fieldChange = -1; %detects field change, set to impossible field number
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip = 5; %number of files to skip
isgood = false(size(datafiles)); %logical array
%loop through all images
for ifile=1:size(datafiles,1)
    
    fprintf('On file %d of %d.\n',ifile,size(datafiles,1));
    
    %load individual mat file for image
    load(sprintf('%s%s',datadir,datafiles(ifile).name));
    
    %skip a # of files at the start
    if( (data.header.fieldnum == fieldChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
        data.header.bad = 1;
        isgood(ifile) = 0; %set logical array
        fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
    elseif( data.header.fieldnum ~= fieldChange )
        fieldChange = data.header.fieldnum;
        fieldChange_fileCntr = 1; %reset
        isgood(ifile) = 0; %set logical array
    else
        data.header.bad = 0;
        isgood(ifile) = 1; %set logical array
    end
    
    im_num(ifile) = ifile;
    exptime(ifile) = data.header.exptime;
    mean_dn_s(ifile) = data.stats.maskmean;
    mean_cal(ifile) = data.stats.calmean;
    mean_dn(ifile) = mean_dn_s(ifile)*exptime(ifile);
    fieldnum(ifile) = data.header.fieldnum;
    header_bias(ifile) = data.ref.biaslevl;
    calc_bias(ifile) = mean(data.ref.line);
    bias_diff(ifile) = header_bias(ifile)-calc_bias(ifile);
    if strcmp(data.ref.biasmthd,'mean of dark column data') == 1
        bias_type(ifile) = 1;   
    elseif strcmp(data.ref.biasmthd,'median of dark column data') == 1
        bias_type(ifile) = 2; 
    else
        fprintf('Absolute panic: bias method not recognized!')
    end
    
end

for i = 1:length(fields)
    j_field = (fieldnum == fields(i)) & isgood;
    field_mean_dn_s(i) = mean(mean_dn_s(j_field));
    field_mean_dn(i) = mean(mean_dn(j_field));
    field_mean_cal(i) = mean(mean_cal(j_field));
end

fprintf('done');