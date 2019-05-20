function cutoff_loop(phase,path)
% load array of file names and JD date
load(sprintf('%smatrix/%s_selected_JD_name.mat',path,phase));
%load(sprintf('%s_selected_JD_name.mat',phase));

% pass filenames into an array for this script
filename = file_array;


median_array = [];
stdev_array=[];

[Nrows, Ncols] = size(filename);

for i=1:Ncols
    if i>(Ncols/4.0) && i < (Ncols/4.0 + 1)
        sprintf('Progress: 25 percent of data processed. Calculating...')
    end
    
    if i>(Ncols/3.0) && i < (Ncols/3.0 + 1)
        sprintf('Progress: 33 percent of data processed. Keep going...')
    end
    
    if i>(Ncols/2.0) && i < (Ncols/2.0 + 1)
        sprintf('Progress: 50 percent of data processed. Looking good...')
    end
    
    if i>(Ncols*0.75) && i < (Ncols*0.75 + 1)
        sprintf('Progress: 75 percent of data processed. Getting there...')
    end
    
    if i>=(Ncols*0.9) && i <= (Ncols*0.9 + 1)
        sprintf('Progress: 75% of data processed. The end is near...')
    end
    
    file_to_read = sprintf('%s%s/selected_data/%s',path,phase,filename(i).name);
    %file_to_read = sprintf('%s/%s',phase,filename(i).name);
    data_raw = single(fitsread(file_to_read));
    data = reshape(data_raw, [numel(data_raw),1]);
    
%     h1 = figure;
%     set(h1,'Visible', 'off');
%     clf;
%     hist(data);
%     set(gca,'YScale','log');
%     ylim([0,100000]);
%     
%     image_old = sprintf('../../data/NH/%s/hist/old_%s',phase,filename(i).name);
%     print(h1,image_old, '-dpng');
%     close(h1);
    
    
    new_data=[];
 
    
    
    
    for k=1:3
       
        standard_dev = std(data);
        median_data = median(data);
        
        for x=1:256
            for y=1:256
            
            data_element = single(data_raw(x,y));
            
                if data_element > (median_data-2.0*standard_dev) && data_element < (median_data+2.0*standard_dev)
                     new_data = [new_data data_element];
                end
        
            end
        end
        
       
%         h2 = figure;
%         set(h2,'Visible', 'off');
%         clf;
%         hist(data);
%         set(gca,'YScale','log');
%         ylim([0,100000]);
%     
%         image_new = sprintf('../../data/NH/%s/hist/%s',phase,filename(i).name);
%         print(h2,image_new, '-dpng');
%         
%         close(h2);
%         
        data = new_data;
        
        size(data);
        
        new_data = 0;
        
%         pause(1);
    end
    
    size(data);
    
    median_array = [median_array median_data];
    stdev_array = [stdev_array standard_dev];
    

end

save(sprintf('%smatrix/%s_med_dev.mat',path,phase), 'date_array', 'file_array', 'median_array', 'stdev_array')
close all
return


