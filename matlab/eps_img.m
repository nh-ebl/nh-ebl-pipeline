function eps_img(phase,min,max,path)

% load date matrix
load(sprintf('%smatrix/%s_time_selected.mat',path, phase));


% load selected files
data = dir(sprintf('%s%s/selected_data/*fit',path,phase));
ndata=numel(data);

% pass the .fit name to filename
filename = sprintf('%s%s/selected_data/*fit',path,phase);
files = dir(filename);
    
for idate= 1:ndata % loop over all selected files
    date = char(date_array_mat(idate));

    file = files(idate).name;
    
    data1=fitsread(sprintf('%s%s/selected_data/%s',path,phase,file));
    
    h = figure(1);
    clf;
    imagesc(data1);
    set(h,'visible','off');
    colormap(flipud(gray));
    colorbar; % display the colorbar on the images
    caxis([min,max]);
    caxis manual;
    
    %add title to the images, which include the date the images were taken
    title(sprintf('%s',date));
    %save as .eps files
    imagename = sprintf('%s%s/eps_images/%s',path,phase,file);
    print(h, imagename, '-deps');
    
    pause(1);
end

return

