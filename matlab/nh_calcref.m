function nh_calcref()

close all

% flag to load previous params, if 'load' loads params, otherwise uses
% whatever string is provided to set params
flag_load_previous = 'new';

if strcmp(flag_load_previous,'load') == 1
    load('run_params.mat','params')
else
    params = struct;
    params.data_type = flag_load_previous; % set the flag's string to be the data_type
    params.err_mags = 0; % turn off
    params.err_str = '';
end
% set paths here for desired data set
if strcmp(params.data_type,'ghost') == 1
    paths = get_paths_old_ghosts();
    ghost = 1;
    old = 0;
    new = 0;
elseif strcmp(params.data_type,'old') == 1
    paths = get_paths_old();
    old = 1;
    ghost = 0;
    new = 0;
elseif strcmp(params.data_type,'new') == 1
    paths = get_paths_new();
    new = 1;
    old = 0;
    ghost = 0;
end

datafiles = dir(sprintf('%s*.mat',paths.datadir));

nfiles = size(datafiles,1);

isgood = zeros(nfiles,1);
mydate = zeros(nfiles,1);
myfieldnum = zeros(nfiles,1);
mytarget = cellstr('');
mysun = zeros(nfiles,1);
mydt = zeros(nfiles,1);
myelong = zeros(nfiles,1);
myexp = zeros(nfiles,1);
myref = zeros(nfiles,1);
myrefold = zeros(nfiles,1);
mynewcorr = zeros(nfiles,1);
myeng = zeros(nfiles,1);
myohm = zeros(nfiles,1);
myav = zeros(nfiles,1);
mysig = zeros(nfiles,1);
myadjmean = zeros(nfiles,1);
myisl = zeros(nfiles,1);
mydgl = zeros(nfiles,1);
mycrr = zeros(nfiles,1);

% set good files based on data set
if new == 1
    goodfiles = [5,6,7,8];
elseif old == 1
    goodfiles = [3,5,6,7];
elseif ghost == 1
    goodfiles = [5,9];
end

for ifile=1:nfiles
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    if params.err_mags == 1
    % Choose mask dir based on err_flag
        datastruct = data.(params.err_str);
    else
        datastruct = data;
    end
    
    mydate(ifile) = data.header.date_jd;
    
    myfieldnum(ifile) = data.header.fieldnum;
    
    mytarget{ifile} = data.header.target_name;
    
    % New reference value - sigma clipped mean
    myref(ifile) = nh_sigclip(datastruct.ref.line);
    % Old reference value - median
    myrefold(ifile) = median(datastruct.ref.line_old);
    % Calculate new reference correction - same as fit that adjusts mean, but
    % applied directly to data
    % corr is later subtracted from data, so this becomes med(ref) is added back,
    % sig_mean(ref) is subtracted, and offset is added
    % Offset determined in nh_dark_analysis_combined from robust fit
    mynewcorr(ifile) = myref(ifile) - myrefold(ifile) - 0.357; 
    % Difference between sigma clipped mean and median
    mycrr(ifile) = datastruct.ref.bias;
    
    % Mean of masked image
    mysig(ifile) = datastruct.stats.maskmean.*data.header.exptime;
    myadjmean(ifile) = mean(datastruct.data(~datastruct.mask.mask) - (mynewcorr(ifile)));
    
    mygood = myfieldnum(ifile) == goodfiles;
    
    if sum(mygood) == 1
        % If struct field data.header.bad exists, check if file is good or
        % bad - only mark as good if not bad
        if isfield(data.header,'bad')
            % If file is not bad, mark as good
            if data.header.bad == 0
                isgood(ifile) = 1;
                
                mystring = sprintf('%d, %s, %f, %f',myfieldnum(ifile),...
                    mytarget{ifile},mysig(ifile),mycrr(ifile));
                
                disp(mystring);
            end
            % If struct field does not exist, files have not been marked good
            % or bad - assume all good
        else
            isgood(ifile) = 1;
            
            mystring = sprintf('%d, %s, %f, %f',myfieldnum(ifile),...
                mytarget{ifile},mysig(ifile),mycrr(ifile));
            
            disp(mystring);
        end
    end
    
end

figure(1)
isgood = logical(isgood);
og = scatter(mycrr(isgood),mysig(isgood),'r');
xlabel('Reference Pixel Bias [DN]');
ylabel('Most Probable Masked Value [DN]');
hold on;
[rho,pval] = corr(mycrr(isgood),mysig(isgood));

[linfit,fiterr] = polyfit(mycrr(isgood),mysig(isgood),1);
ste = sqrt(diag(inv(fiterr.R)*inv(fiterr.R')).*fiterr.normr.^2./fiterr.df);
xcrr=linspace(min(mycrr(isgood)),max(mycrr(isgood)));
sigfit=(linfit(1).*xcrr + linfit(2));
fit = plot(xcrr,sigfit);

% Original fit - makes mean 0
% mycorr = linfit(1).*mycrr + linfit(2);
% New fit - makes mean be previous y-intercept
mycorr = linfit(1).*mycrr + linfit(2) - linfit(2);

[rho,pval] = corr(mycrr(isgood),mysig(isgood)-mycorr(isgood));

cor = scatter(mycrr(isgood),mysig(isgood)-mycorr(isgood),'b');
% cor = scatter(mycrr(isgood),mycorr(isgood),'b');
% title('calcref');
xlabel('Reference Pixel Bias [DN]');
% ylabel('Most Probable Masked Value [DN]');
ylabel('Mean of Masked Image [DN]');
% legend([og,cor,fit],{'Original','Original - Corrected','Linear Fit'},'Location','northeast');
legend([og,cor,fit],{'Original','Corrected','Linear Fit'},'Location','northeast');
correction.date = mydate;
correction.corr = mycorr;

% Plot comparison between old ref (med(ref)) with old masked mean (based on
% med(ref) subtraction), new ref (sigma-mean(ref)) with new masked mean
% (based on sigma-mean(ref) subtraction, and new ref with fit-corrected old
% masked mean
figure(2)
% Calculate fit for desired pairing
[reffit,reffiterr] = polyfit(myref(isgood),myadjmean(isgood),1);
xref=linspace(min(myref(isgood)),max(myref(isgood)));
meanfit=(reffit(1).*xref + reffit(2));
og = scatter(myrefold(isgood),mysig(isgood),'r');
hold on;
ne = scatter(myref(isgood),myadjmean(isgood),'b');
fix = scatter(myref(isgood),mysig(isgood)-mycorr(isgood),'g');
fit2 = plot(xref,meanfit);
xlabel('Reference Med/Mean [DN]');
ylabel('Mean of Masked Image [DN]');
legend([og,ne,fix],{'Med Ref','Sig Ref with Offset','Adj Med Ref'},'Location','northeast');
title(sprintf('Fit: y = %.2fx + %.2f',reffit(1),reffit(2))); 

if new == 1
    save('lookup/nh_refcorr_new.mat','correction');
    % Save correction to each data file
    for ifile=1:nfiles
        load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
        
        if params.err_mags == 1
            datastruct = data.(params.err_str);
        else
            datastruct = data;
        end
        
        % Save new version based on most probable value of masked image
        %         data.refcorr.mostprob_date = mydate;
        %         data.refcorr.mostprob_corr = mycorr;
        % Save original (based on masked image mean)
        datastruct.refcorr.date = mydate;
        datastruct.refcorr.corr = mycorr;
        
        if params.err_mags == 1
            data.(params.err_str) = datastruct;
        else
            data = datastruct;
        end
        
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
    end
elseif old == 1
    save('lookup/nh_refcorr_old.mat','correction');
elseif ghost == 1
    save('lookup/nh_refcorr_ghost.mat','correction');
end

end