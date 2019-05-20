%This a function to create png images from fits images
%It includes the option to input min and max color scaling values or to use
%automatic scaling

%Symons, 2018

function png_img_cust(phase,min,max,path)

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
%     caxis([min,max]);
%     caxis manual;
    caxis('auto');
    
    %add title to the images, which include the date the images were taken
    title(sprintf('%s',file));
    %save as .eps files
    imagenamefits = sprintf('%s%s/check_files/%s',path,phase,file);
    imagenamepng = strrep(imagenamefits,'fit','png')
    print(h,imagenamepng, '-dpng');
    
    pause(1);
end

return

