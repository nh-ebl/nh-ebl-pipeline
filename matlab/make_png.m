clear all 
close all

path = '/data/symons/nh_data/'
phase = 'pluto'

% load date matrix
load(sprintf('%smatrix/%s_time_selected.mat',path, phase));


% load selected files
data = dir(sprintf('%s%s/check_files/*fit',path,phase));
ndata=numel(data);

% pass the .fit name to filename
filename = sprintf('%s%s/check_files/*fit',path,phase);
files = dir(filename);
    
for idate= 1:ndata % loop over all selected files
    date = char(date_array_mat(idate));

    file = files(idate).name;
    
    data1=fitsread(sprintf('%s%s/check_files/%s',path,phase,file));
    
    h = figure(1);
    clf;
    imagesc((data1));
    set(h,'visible','off');
    colormap(flipud(gray));
    colorbar; % display the colorbar on the images
    caxis([-5,50]);
    caxis manual;
    
    %add title to the images, which include the date the images were taken
    title(sprintf('%s',date));
    %save as .eps files
    ext = '.png'
    imagename = sprintf('%s%s/check_files/%s%s',path,phase,file,ext);
    print(h,imagename, '-dpng');
    
    pause(1);
end


