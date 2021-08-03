 function [data, evenMinusOdd, evenMinusOddFixed, oddMinusEvenref, oddMinusEvenrefFixed] = nh_jail_bars(data,paths,flag_method)
%% retrieve original data.data from fits file
% read in the .fits data for the calibrated image file
imagename = ['regist_',data.header.rawfile];
image_i = fitsread(sprintf('%s%s',paths.imagedir,imagename));

% append the image data to the data structure
data.data = image_i;

%save old data.data for reference - can revert to this if needed
data.image.original_data = data.data;

% read in the .fits data for the raw (engineering) image file - includes reference column
ref_i = fitsread(sprintf('%s',data.ref.file));
data.ref.line_old = ref_i(:,257);

% If not using new method, skip doing jail bar correction
% if (strcmp(flag_method, 'old_corr') == 1 || strcmp(flag_method,'new') == 1 || strcmp(flag_method,'old') == 1)
if (strcmp(flag_method,'new') == 1)

    evenMinusOdd = 1;
    evenMinusOddFixed = 1;

    %% checking LORRI images for jail bars and fixing if necessary

    imageNaN = data.data.*~data.mask.onemask; %get image
    imageNaN(data.mask.onemask) = nan; %NaN masked pixels (so mean doesn't involve them)

    maskvect = ones([256,1],'logical');
    refmask = horzcat(data.mask.onemask, maskvect);
    refNaN = ref_i.*~refmask; %get image
    refNaN(refmask) = nan; %NaN masked pixels (so mean doesn't involve them)

    oddvar = imageNaN(:,1:2:end); %get odd columns
    evenvar = imageNaN(:,2:2:end); %get even columns

    oddvarref = refNaN(:,3:2:end); %get odd columns
    evenvarref = refNaN(:,2:2:end); %get even columns

    %subtract columns, 1st 'nanmean of columns' then 2nd nanmean the 'nanmean of columns'
    evenMinusOdd = nanmean(nanmean(evenvar-oddvar));
    oddMinusEvenref = nanmean(nanmean(oddvarref-evenvarref));

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


    if(oddMinusEvenref > 0)
        %odd columns have the jail bars, subtract them off
        ref_i(:,1:2:end) = ref_i(:,1:2:end) - abs(oddMinusEvenref);
    else
        %even columns have the jail bars, subtract them off
        ref_i(:,2:2:end) = ref_i(:,2:2:end) - abs(oddMinusEvenref);
    end

    %save jail bar adjustment for reference
    data.header.jailbar_adj = abs(evenMinusOdd);
    data.header.jailbar_adj_ref = abs(oddMinusEvenref);

    %redo it all for testing
    imageNaN = data.data.*~data.mask.onemask;
    imageNaN(data.mask.onemask) = nan;
    oddvar = imageNaN(:,1:2:end);
    evenvar = imageNaN(:,2:2:end);
    evenMinusOddFixed = nanmean(nanmean(evenvar-oddvar));

    refNaN = ref_i.*~refmask; %get image
    refNaN(refmask) = nan; %NaN masked pixels (so mean doesn't involve them)
    oddvarref = refNaN(:,3:2:end); %get odd columns
    evenvarref = refNaN(:,2:2:end); %get even columns
    oddMinusEvenrefFixed = nanmean(nanmean(oddvarref-evenvarref));
else
    evenMinusOdd = 0;
    evenMinusOddFixed = 0;
    oddMinusEvenref = 0;
    oddMinusEvenrefFixed = 0;
end

% recalculate mask mean and std with new data.data
datmean = mean(data.data(~data.mask.mask)./ data.astrom.exptime);
datstd = std(data.data(~data.mask.mask)./ data.astrom.exptime);

% and append it to the data structure
data.stats.maskmean = datmean;
data.stats.maskstd = datstd;
data.stats.maskerr = datstd ./ sqrt(256.^2 - sum(data.mask.onemask(:)));

% recalculate reference image and column info and save
data.ref.eng = ref_i;
data.ref.line = ref_i(:,257);
whpl = ~isnan(data.ref.line);
data.ref.mean = mean(ref_i(whpl,257));
data.ref.std = std(ref_i(whpl,257));

temp_ref = data.ref.eng(:,1:256);

data.ref.engmean = mean(temp_ref(~data.mask.onemask));
data.ref.bias = nh_sigclip(data.ref.line) - median(data.ref.line_old);

end

