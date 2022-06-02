% This script independently analyzes the startup masked image mean of each
% sequence of observations to look for a positive bias due to camera
% activation

% Can examine any or all data sets

clear all
close all
clc

%Import paths for data location
paths = get_paths_old();
npaths = get_paths_new();
lpaths = get_paths_lauer();
wpaths = get_paths_newest();

%old light data
lightfiles = dir(sprintf('%s*.mat',paths.datadir));
%new light data
nlightfiles = dir(sprintf('%s*.mat',npaths.datadir));
%lauer light data
llightfiles = dir(sprintf('%s*.mat',lpaths.datadir));
%newest light data
wlightfiles = dir(sprintf('%s*.mat',wpaths.datadir));

%Determine which light files correspond to 'good' fields
isgoodold = zeros(numel(lightfiles),1);
isgoodnew = zeros(numel(nlightfiles),1);
isgoodlauer = zeros(numel(llightfiles),1);
isgoodnewest = zeros(numel(wlightfiles),1);

%These field numbers were previously determined by going through all the
%files
oldgoodfields = [3,5,6,7];
% oldgoodfields = []; % Set good fields to none to skip data set
newgoodfields = [5,6,7,8];
% newgoodfields = [];
lauergoodfields = [1,2,3,4,5,6,7];
% lauergoodfields = [];
lauer_exlude_enable = false; %enables skipping of new sequences
newestgoodfields = [2,4,5,6,7,12,15,16,17,19,20,22,23];
% newestgoodfields = [];
newest_exlude_enable = false; %enables skipping of new sequences

%Check for old light files
parfor ifile=1:numel(lightfiles)
    datatemp = load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == oldgoodfields)
        if isfield(data.header,'bad')
            if data.header.bad == 0
                isgoodold(ifile) = 1;
            end
        else
            isgoodold(ifile) = 1;
        end
    end
end

%Check for new light files
parfor ifile=1:numel(nlightfiles)
    datatemp = load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == newgoodfields)
        if isfield(data.header,'bad')
            if data.header.bad == 0
                isgoodnew(ifile) = 1;
            end
        else
            isgoodnew(ifile) = 1;
        end
    end
end

%Check for lauer light files
parfor ifile=1:numel(llightfiles)
    datatemp = load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == lauergoodfields)
        if isfield(data.header,'bad')
            if data.header.bad == 0
                isgoodlauer(ifile) = 1;
            end
        else
            isgoodlauer(ifile) = 1;
        end
    end
end

%Check for newest light files
parfor ifile=1:numel(wlightfiles)
    datatemp = load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
    data = datatemp.data; %allows parallel to work
    if sum(data.header.fieldnum == newestgoodfields)
        if isfield(data.header,'bad')
            if data.header.bad == 0
                isgoodnewest(ifile) = 1;
            end
        else
            isgoodnewest(ifile) = 1;
        end
    end
end

%Number of old and new light files corresponding to good fields
numoldlightfiles = sum(isgoodold);
numnewlightfiles = sum(isgoodnew);
numlauerlightfiles = sum(isgoodlauer);
numnewestlightfiles = sum(isgoodnewest);

%Preallocate
im_mean_old = zeros((numoldlightfiles),3);
im_mean_old_names = cell(26,1);
im_mean_new = zeros((numnewlightfiles),3);
im_mean_new_names = cell(343,1);
im_mean_lauer = zeros((numlauerlightfiles),3);
im_mean_lauer_names = cell(320,1);
im_mean_newest = zeros((numnewestlightfiles),3);
im_mean_newest_names = cell(1106,1);

% Reset sequence info
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file keep
fieldChange_fileSkip_time = 300; %s, time to keep at start of sequence
FLG_firstRun = true; %set first run flag
reqIDCounter = 1; %counter for number of unique reqIDs

%For old data files
fprintf('Loading old data \n');
for ifile=1:numel(lightfiles)
    %If file is for a good field, load and save values
    if isgoodold(ifile) == 1
        %Load data files
        load(sprintf('%s%s',paths.datadir,lightfiles(ifile).name));
        if( FLG_firstRun )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to keep dynamically
            FLG_firstRun = false; %turn off
        end
        % Save means for a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            im_mean_old(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_old(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_old_names{ifile,1} = reqIDChange; %record reqID
            im_mean_old(ifile,3) = data.astrom.exptime; %set exptime
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            im_mean_old(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_old_names{ifile,1} = reqIDChange; %record reqID
            if ifile > 1
                reqIDCounter = reqIDCounter + 1; % increment
            end
            im_mean_old(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_old(ifile,3) = data.astrom.exptime; %set exptime
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to keep
        end
    end
end

% Reset sequence info
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file keep
fieldChange_fileSkip_time = 300; %s, time to keep at start of sequence
FLG_firstRun = true; %set first run flag
reqIDCounter = 1; %counter for number of unique reqIDs

%For new data files
fprintf('Loading new data \n');
for ifile=1:numel(nlightfiles)
    %If file is for a good field, load and save values
    if isgoodnew(ifile) == 1
        %Load data files
        load(sprintf('%s%s',npaths.datadir,nlightfiles(ifile).name));
        if( FLG_firstRun )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to keep dynamically
            FLG_firstRun = false; %turn off
        end
        % Save means for a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            im_mean_new(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_new(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_new_names{ifile,1} = reqIDChange; %record reqID
            im_mean_new(ifile,3) = data.astrom.exptime; %set exptime
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            im_mean_new(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_new_names{ifile,1} = reqIDChange; %record reqID
            if ifile > 1
                reqIDCounter = reqIDCounter + 1; % increment
            end
            im_mean_new(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_new(ifile,3) = data.astrom.exptime; %set exptime
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to keep
        end
    end
end

% Reset sequence info
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file keep
fieldChange_fileSkip_time = 300; %s, time to keep at start of sequence
FLG_firstRun = true; %set first run flag
reqIDCounter = 1; %counter for number of unique reqIDs

%For lauer data files
fprintf('Loading Lauer data \n');
jfile = 1;
for ifile=1:numel(llightfiles)
    %If file is for a good field, load and save values
    if isgoodlauer(ifile) == 1
        %Load data files
        load(sprintf('%s%s',lpaths.datadir,llightfiles(ifile).name));
        if( FLG_firstRun )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to keep dynamically
            FLG_firstRun = false; %turn off
        end
        % Save means for a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            im_mean_lauer(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_lauer(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_lauer_names{ifile,1} = reqIDChange; %record reqID
            im_mean_lauer(ifile,3) = data.astrom.exptime; %set exptime
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            im_mean_lauer(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_lauer_names{ifile,1} = reqIDChange; %record reqID
            if ifile > 1
                reqIDCounter = reqIDCounter + 1; % increment
            end
            im_mean_lauer(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_lauer(ifile,3) = data.astrom.exptime; %set exptime
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to keep
        end
    end
end

% Reset sequence info
reqIDChange = ''; %detects reqID change
fieldChange_fileCntr = 1; %counter for file keep
fieldChange_fileSkip_time = 300; %s, time to keep at start of sequence
FLG_firstRun = true; %set first run flag
reqIDCounter = 1; %counter for number of unique reqIDs

%For newest data files
fprintf('Loading newest data \n');
jfile = 1;
for ifile=1:numel(wlightfiles)
    %If file is for a good field, load and save values
    if isgoodnewest(ifile) == 1
        %Load data files
        load(sprintf('%s%s',wpaths.datadir,wlightfiles(ifile).name));
        if( FLG_firstRun )
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %get how many files to keep dynamically
            FLG_firstRun = false; %turn off
        end
        % Save means for a # of files at the start
        if( strcmp(data.astrom.reqid, reqIDChange) && (fieldChange_fileCntr < fieldChange_fileSkip) )
            fieldChange_fileCntr = fieldChange_fileCntr + 1; %increment
            im_mean_newest(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_newest(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_newest_names{ifile,1} = reqIDChange; %record reqID
            im_mean_newest(ifile,3) = data.astrom.exptime; %set exptime
        elseif( ~strcmp(data.astrom.reqid, reqIDChange) )
            reqIDChange = data.astrom.reqid; %record reqID
            fieldChange_fileCntr = 1; %reset
            im_mean_newest(ifile,1) = data.stats.maskmean; %record image mean
            im_mean_newest_names{ifile,1} = reqIDChange; %record reqID
            if ifile > 1
                reqIDCounter = reqIDCounter + 1; % increment
            end
            im_mean_newest(ifile,2) = reqIDCounter; %set current number of reqID
            im_mean_newest(ifile,3) = data.astrom.exptime; %set exptime
            fieldChange_fileSkip = round(fieldChange_fileSkip_time/data.header.exptime); %recalc how many fields to keep
        end
    end
end

% Remove zero entries
old_means = im_mean_old(abs(im_mean_old(:,1)) > 0,:);
new_means = im_mean_new(abs(im_mean_new(:,1)) > 0,:);
lauer_means = im_mean_lauer(abs(im_mean_lauer(:,1) )> 0,:);
newest_means = im_mean_newest(abs(im_mean_newest(:,1) )> 0,:);

%old
curr_means = old_means; %save rewriting code
curr_title = 'Old Images'; %save rewriting code
[curr_seqs, ~, curr_seqsIndx] = unique(curr_means(:,2)); %get number of sequencies
curr_lines = []; %holder for line handles
curr_legend = cell(1,length(curr_seqs)); %holder for legend text
figure; clf
hold on;
for i = 1:length(curr_seqs)
    kj = i == curr_seqsIndx; %get current sequence location
    curr_lines(end+1) = plot(cumsum(curr_means(kj,3)), curr_means(kj,1),'-o'); %plot and get handle
    curr_legend{1,i} = num2str(curr_seqs(i)); %record legend text
end
hLeg = legend(curr_lines,curr_legend); %make a legend
% set(hLeg,'visible','off'); %turn off legend, just needed it for names
ylim([0 0.11])
xlabel('Time Since Sequence Start (s)');
ylabel('Mean Sky Level [DN/s]');
title(curr_title);

%new
curr_means = new_means; %save rewriting code
curr_title = 'New Images'; %save rewriting code
[curr_seqs, ~, curr_seqsIndx] = unique(curr_means(:,2)); %get number of sequencies
curr_lines = []; %holder for line handles
curr_legend = cell(1,length(curr_seqs)); %holder for legend text
figure; clf
hold on;
for i = 1:length(curr_seqs)
    kj = i == curr_seqsIndx; %get current sequence location
    curr_lines(end+1) = plot(cumsum(curr_means(kj,3)), curr_means(kj,1),'-o'); %plot and get handle
    curr_legend{1,i} = num2str(curr_seqs(i)); %record legend text
end
hLeg = legend(curr_lines,curr_legend); %make a legend
% set(hLeg,'visible','off'); %turn off legend, just needed it for names
ylim([0 0.11])
xlabel('Time Since Sequence Start (s)');
ylabel('Mean Sky Level [DN/s]');
title(curr_title);

%lauer
curr_means = lauer_means; %save rewriting code
curr_title = 'Lauer Images'; %save rewriting code
[curr_seqs, ~, curr_seqsIndx] = unique(curr_means(:,2)); %get number of sequencies
curr_lines = []; %holder for line handles
curr_legend = cell(1,length(curr_seqs)); %holder for legend text
hFig = figure; clf
datacursormode on;
dcm = datacursormode(hFig);
set(dcm,'UpdateFcn',@dataTipCustom) %for legend off
hold on;
for i = 1:length(curr_seqs)
    kj = i == curr_seqsIndx; %get current sequence location
    curr_lines(end+1) = plot(cumsum(curr_means(kj,3)), curr_means(kj,1),'-o'); %plot and get handle
    curr_legend{1,i} = num2str(curr_seqs(i)); %record legend text
end
hLeg = legend(curr_lines,curr_legend); %make a legend
set(hLeg,'visible','off'); %turn off legend, just needed it for names
ylim([0 0.11])
xlabel('Time Since Sequence Start (s)');
ylabel('Mean Sky Level [DN/s]');
title(curr_title);

%newest
curr_means = newest_means; %save rewriting code
curr_title = 'Newest Images'; %save rewriting code
[curr_seqs, ~, curr_seqsIndx] = unique(curr_means(:,2)); %get number of sequencies
curr_lines = []; %holder for line handles
curr_legend = cell(1,length(curr_seqs)); %holder for legend text
hFig2 = figure; clf
datacursormode on;
dcm2 = datacursormode(hFig2);
set(dcm2,'UpdateFcn',@dataTipCustom) %for legend off
hold on;
for i = 1:length(curr_seqs)
    kj = i == curr_seqsIndx; %get current sequence location
    curr_lines(end+1) = plot(cumsum(curr_means(kj,3)), curr_means(kj,1),'-o'); %plot and get handle
    curr_legend{1,i} = num2str(curr_seqs(i)); %record legend text
end
hLeg = legend(curr_lines,curr_legend); %make a legend
set(hLeg,'visible','off'); %turn off legend, just needed it for names
% ylim([0 0.11])
xlabel('Time Since Sequence Start (s)');
ylabel('Mean Sky Level [DN/s]');
title(curr_title);


fprintf('done');



function output_txt = dataTipCustom(obj,event_obj,str)
    %for legends that are too big
    pos = get(event_obj, 'Position');
    output_txt = {...
        ['X: ', num2str(pos(1),4)]...
        ['Y: ', num2str(pos(2),4)] ...
        ['#:', event_obj.Target.DisplayName]...
    };
end