function nh_calcref()

% set paths here for desired data set
% paths = get_paths_old_ghosts();
paths = get_paths_old();
% paths = get_paths_new();

% based on paths, decide if we're analysing old data, new data, or ghost
% data
if strcmp(paths.datadir,'/data/symons/NH_old_data/mat/ghosts/') == 1
    ghost = 1;
    old = 0;
    new = 0;
elseif strcmp(paths.datadir,'/data/symons/NH_old_data/mat/good/') == 1
    old = 1;
    ghost = 0;
    new = 0;
elseif strcmp(paths.datadir,'/data/symons/nh_data/mat/') == 1
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
myeng = zeros(nfiles,1);
myohm = zeros(nfiles,1);
myav = zeros(nfiles,1);
mysig = zeros(nfiles,1);
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
    
    mydate(ifile) = data.header.date_jd;
    
    myfieldnum(ifile) = data.header.fieldnum;
    
    mytarget{ifile} = data.header.target_name;
    
    % Mean of masked image
    mysig(ifile) = data.stats.maskmean.*data.header.exptime;
    % Most probable masked value
    %     mysig(ifile) = data.stats.mostprob.*data.header.exptime;
    
    mycrr(ifile) = data.ref.bias;
    
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

mycorr = linfit(1).*mycrr + linfit(2);

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

if new == 1
    save('lookup/nh_refcorr_new.mat','correction');
    % Save correction to each data file
    for ifile=1:nfiles
        load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
        % Save new version based on most probable value of masked image
        %         data.refcorr.mostprob_date = mydate;
        %         data.refcorr.mostprob_corr = mycorr;
        % Save original (based on masked image mean)
        data.refcorr.date = mydate;
        data.refcorr.corr = mycorr;
        save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
    end
elseif old == 1
    save('lookup/nh_refcorr_old.mat','correction');
elseif ghost == 1
    save('lookup/nh_refcorr_ghost.mat','correction');
end

end