% This script tests various values of N used for an N*sigma clip in the
% creation of the clip mask, performed in nh_makemask. For all available
% data files, plots of mu after masking versus every value of N are made
% (or alternatively overplotted on a single plot). The derivative of mu
% versus N can also be calculated and plotted, or an average derivative for
% all data files can be plotted. This is used to determine the ideal N
% value for the clip mask. 

% Symons 2019

clear all
close all

% Retrieve paths for new data
paths = get_paths_new();

% Create list of .mat files
datafiles = dir(sprintf('%s*.mat',paths.datadir));

% Define N for N*sigma used in clip mask
nsig = linspace(0.1,5,50);

% Build colormap
cmap = colormap(parula(size(datafiles,1))); %blue to yellow
% cmap = colormap(hsv(size(datafiles,1))); %rainbow, red to red

% Preallocate
grad = zeros(length(nsig),size(datafiles,1));
dydx = zeros(length(nsig),size(datafiles,1));

% For all data files
for ifile=1:size(datafiles,1)
    
    % Preallocate for mu
    checkmu = zeros(length(nsig),1);
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    % Load data file
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    % For all test values of N
    for i=1:length(nsig)
        
        % Run makemask with current N value and retrieve mu
        [data,mu] = nh_makemask(data,paths,nsig(i));
        checkmu(i) = mu;
    end
    
    % Calculate gradient and derivative of mu with respect to N
    grad(:,ifile) = gradient(checkmu)./gradient(nsig');
    dydx(:,ifile) = diff([eps;checkmu])./diff([eps;nsig']);
    
    % Plot mu versus N for all data files
%     h = figure(1);
%     hold on
%     plot(nsig',grad,'color',cmap(ifile,:));
% %     scatter(nsig,checkmu,5,cmap(ifile,:),'filled');
%     xlabel('N sigma for clip mask');
%     ylabel('Mu after mask');
    
end

% Plot average derivative of mu versus N
h = figure(1);
plot(nsig',mean(dydx,2));
xlabel('N sigma for clip mask');
ylabel('Mu after mask');
