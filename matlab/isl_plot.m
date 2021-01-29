clear all
close all

mag = 16.5;

filename = strcat(num2str(mag),'_mask.mat');
load(filename);
mask_list = star_list;

filename = strcat(num2str(mag),'_gaia.mat');
load(filename);
gaia_list = star_list;

filename = strcat(num2str(mag),'_tri_only.mat');
load(filename);
tri_only_list = star_list;

filename = strcat(num2str(mag),'_wing.mat');
load(filename);
wing_list = star_list;

load('20_tri.mat');
tri_list = star_list;

h = figure;
clf;
hold on
scatter(mask_list(:,1),mask_list(:,2),70,'filled');
scatter(wing_list(:,1),wing_list(:,2),'filled');
scatter(tri_only_list(:,1),tri_only_list(:,2),70,'filled');
scatter(tri_list(:,1),tri_list(:,2),50,'filled');
scatter(gaia_list(:,1),gaia_list(:,2),30,'filled');
set(gca,'YScale','log');
legend('Masked < m','Wing < m','Trilegal Only > m','Trilegal > 20','m < Gaia < 20');
title(sprintf('m = %1.1f',mag));
xlabel('Mag of Star');
ylabel('Flux Contribution [nW m^-2 sr^-1]');
% ylim( [0.5, 100000] ); %force 0
% xlim([-50,50]);
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.histdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

magBinsStepSize = 1; %sets the step size for the mag bins
magBinsCat = horzcat((min(floor(mask_list(:,1))):magBinsStepSize:max(ceil(mask_list(:,1)))-magBinsStepSize)',(min(floor(mask_list(:,1)))+magBinsStepSize:magBinsStepSize:max(ceil(mask_list(:,1))))'); %get the mag bins spaced 1 apart
magBinsSizeCat = size(magBinsCat); %get size of mag bins since can't compound calls in matlab
magBinsTri = horzcat((min(floor(tri_only_list(:,1))):magBinsStepSize:max(ceil(tri_only_list(:,1)))-magBinsStepSize)',(min(floor(tri_only_list(:,1)))+magBinsStepSize:magBinsStepSize:max(ceil(tri_only_list(:,1))))'); %get the mag bins spaced 1 apart
magBinsSizeTri = size(magBinsTri); %get size of mag bins since can't compound calls in matlab

for irate = 1:magBinsSizeCat(1)
    k = (mask_list(:,1) >= magBinsCat(irate,1)) & (mask_list(:,1) < magBinsCat(irate,2)); %get where magnitudes are within the bin
    curr_mymag_cat = mask_list(k,1); %get only values related to the mag range
    curr_myflux_cat = mask_list(k,2); %get only values related to the mag range
    starNumcat = length(curr_mymag_cat); %num of stars from data
    catnum(irate) = starNumcat; %record the number of stars in bin
    catsum(irate) = sum(curr_myflux_cat); %record the total flux of stars in bin
end

for irate = 1:magBinsSizeTri(1)
    k = (tri_only_list(:,1) >= magBinsTri(irate,1)) & (tri_only_list(:,1) < magBinsTri(irate,2)); %get where magnitudes are within the bin
    curr_mymag_tri = tri_only_list(k,1); %get only values related to the mag range
    curr_myflux_tri = tri_only_list(k,2); %get only values related to the mag range
    starNumtri = length(curr_mymag_tri); %num of stars from data
    trinum(irate) = starNumtri; %record the number of stars in bin
    trisum(irate) = sum(curr_myflux_tri); %record the total flux of stars in bin
end

h = figure;
clf;
hold on
scatter(magBinsCat(:,1)+0.5,catsum,70,'filled');
scatter(magBinsTri(2:end,1)+0.5,trisum(:,2:end),70,'filled');
%     scatter(wing_list(:,1),wing_list(:,2),'filled');
%     scatter(tri_only_list(:,1),tri_only_list(:,2),70,'filled');
%     scatter(tri_list(:,1),tri_list(:,2),50,'filled');
%     scatter(gaia_list(:,1),gaia_list(:,2),30,'filled');
set(gca,'YScale','log');
legend('Masked < m','Trilegal Only > m');
title(sprintf('m = %1.1f',mag));
xlabel('Mag Bin');
ylabel('Surface Brightness in Bin [nW m^-2 sr^-1]');