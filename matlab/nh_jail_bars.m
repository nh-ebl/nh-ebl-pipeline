function [data, evenMinusOdd, evenMinusOddFixed] = nh_jail_bars(data)


%get paths for new data files or old data files
paths = get_paths_new();

%% retrieve original data.data from fits file
% read in the .fits data for the file 
imagename = ['regist_',data.header.rawfile];
image_i = fitsread(sprintf('%s%s',paths.imagedir,imagename));

% append the image data to the data structure
data.data = image_i;

%save old data.data for reference - can revert to this if needed
data.image.original_data = data.data;

evenMinusOdd = 1;
evenMinusOddFixed = 1;

%% checking LORRI images for jail bars and fixing if necessary

imageNaN = data.data.*~data.mask.onemask; %get image
imageNaN(data.mask.onemask) = nan; %NaN masked pixels (so mean doesn't involve them)

oddvar = imageNaN(:,1:2:end); %get odd columns
evenvar = imageNaN(:,2:2:end); %get even columns

%subtract columns, 1st 'nanmean of columns' then 2nd nanmean the 'nanmean of columns'
evenMinusOdd = nanmean(nanmean(evenvar-oddvar)); 

if(evenMinusOdd > 0)
    %even columns have the jail bars, subtract them off
    data.data(:,2:2:end) = data.data(:,2:2:end) - abs(evenMinusOdd);
    %save jail bar type
    data.header.jailbar_type = 2;
else
    %odd columns have the jail bars, subtract them off
    data.data(:,1:2:end) = data.data(:,1:2:end) - abs(evenMinusOdd);
    %save jail bar type
    data.header.jailbar_type = 1;
end

%save jail bar adjustment for reference
data.header.jailbar_adj = abs(evenMinusOdd);

%redo it all for testing
imageNaN = data.data.*~data.mask.onemask;
imageNaN(data.mask.onemask) = nan;
oddvar = imageNaN(:,1:2:end);
evenvar = imageNaN(:,2:2:end);
evenMinusOddFixed = nanmean(nanmean(evenvar-oddvar));

% recalculate mask mean and std with new data.data
datmean = mean(data.data(~data.mask.mask)./ data.astrom.exptime);
datstd = std(data.data(~data.mask.mask)./ data.astrom.exptime);

% and append it to the data structure
data.stats.maskmean = datmean;
data.stats.maskstd = datstd;
data.stats.maskerr = datstd ./ sqrt(256.^2 - sum(data.mask.onemask(:)));

end

