
clear all
close all

paths = get_paths_new();

datafiles = dir(sprintf('%s*.mat',paths.datadir));

nsig = linspace(0.1,5,50);

% Build colormap
cmap = colormap(parula(size(datafiles,1))); %blue to yellow
% cmap = colormap(hsv(size(datafiles,1))); %rainbow, red to red

grad = zeros(length(nsig),size(datafiles,1));
dydx = zeros(length(nsig),size(datafiles,1));

for ifile=1:size(datafiles,1)
    
    checkmu = zeros(length(nsig),1);
    
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    for i=1:length(nsig)
    
        [data,mu] = nh_makemask(data,paths,nsig(i));
        checkmu(i) = mu;
    end
    
    grad(:,ifile) = gradient(checkmu)./gradient(nsig');
    dydx(:,ifile) = diff([eps;checkmu])./diff([eps;nsig']);
    
%     h = figure(1);
%     hold on
%     plot(nsig',grad,'color',cmap(ifile,:));
% %     scatter(nsig,checkmu,5,cmap(ifile,:),'filled');
%     xlabel('N sigma for clip mask');
%     ylabel('Mu after mask');
    
end


h = figure(1);
plot(nsig',mean(dydx,2));
xlabel('N sigma for clip mask');
ylabel('Mu after mask');
