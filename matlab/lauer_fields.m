clear all
close all

%get paths for new data files or old data files
paths = get_paths_lauer();
datadir = paths.datadir;

%load directory of saved mat files for all images
datafiles = dir(sprintf('%s*.mat',datadir));
imagedata = dir(sprintf('%s*.fit',paths.imagedir));

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
field_mean_fieldChangeOnly_dn_s = zeros(size(fields));
field_mean_fieldChangeOnly_dn = zeros(size(fields));
field_mean_fieldChangeOnly_cal = zeros(size(fields));
field_dn_s = zeros(size(fields));
field_dn = zeros(size(fields));
field_cal = zeros(size(fields));


reqIDChange = ''; %detects reqID change
badFileCntr = 0; %counts bad fields
fieldChange = -1; %detects field change, set to impossible field number
fieldChange_fileCntr = 1; %counter for file skip
fieldChange_fileSkip = 5; %number of files to skip
isgood = false(size(datafiles)); %logical array
isgood_fieldChange = false(size(datafiles)); %logical array
%loop through all images
for ifile=1:size(datafiles,1)
    
    fprintf('On file %d of %d.\n',ifile,size(datafiles,1));
    
    %load individual mat file for image
    load(sprintf('%s%s',datadir,datafiles(ifile).name));

    info = fitsinfo([paths.imagedir,imagedata(ifile).name]);
    
    image_i = fitsread(sprintf('%s%s',...
        paths.imagedir,imagedata(ifile).name));
    
    %this is to see what fieldnum isgood looked like compared to reqid
    if( (data.header.fieldnum == fieldChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
        isgood_fieldChange(ifile) = 0; %set logical array
    elseif( data.header.fieldnum ~= fieldChange )
        fieldChange = data.header.fieldnum;
        isgood_fieldChange(ifile) = 0; %set logical array
    else
        isgood_fieldChange(ifile) = 1; %set logical array
    end
    %skip a # of files at the start
    if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
        data.header.bad = 1;
        isgood(ifile) = 0; %set logical array
        fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
        badFileCntr = badFileCntr + 1; %increment
    elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
        %fieldChange = data.header.fieldnum;
        reqIDChange = data.astrom.reqid; %record reqID
        fieldChange_fileCntr = 1; %reset
        data.header.bad = 1;
        isgood(ifile) = 0; %set logical array
        badFileCntr = badFileCntr + 1; %increment
    else
        data.header.bad = 0;
        isgood(ifile) = 1; %set logical array
    end
    
    im_num(ifile) = ifile;
    exptime(ifile) = data.header.exptime;
    mean_dn_s(ifile) = data.stats.maskmean;
    mean_cal(ifile) = data.stats.calmean;
    mean_dn(ifile) = mean_dn_s(ifile)*exptime(ifile);
    disp([data.header.timestamp,' mean: ',num2str(mean_dn(ifile)),' skipped: ',num2str(~isgood(ifile))])
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

disp(['In lauer_fields.m: skipped ',num2str(badFileCntr),' bad files.']) %report bad fields skipped

for i = 1:length(fields)
    j_field = (fieldnum == fields(i)) & isgood;
    field_mean_dn_s(i) = mean(mean_dn_s(j_field));
    field_mean_dn(i) = mean(mean_dn(j_field));
    field_mean_cal(i) = mean(mean_cal(j_field));

    j_field = (fieldnum == fields(i)) & isgood_fieldChange;
    field_mean_fieldChangeOnly_dn_s(i) = mean(mean_dn_s(j_field));
    field_mean_fieldChangeOnly_dn(i) = mean(mean_dn(j_field));
    field_mean_fieldChangeOnly_cal(i) = mean(mean_cal(j_field));
end

fprintf(['field_mean_dn REQID CHANGE:\t\t',num2str(field_mean_dn),'\nfield_mean_dn FIELDNUM CHANGE (old):\t',num2str(field_mean_fieldChangeOnly_dn),...
    '\nDiff:\t\t\t\t\t',num2str(field_mean_dn-field_mean_fieldChangeOnly_dn),'\n'])

fprintf('done\n');