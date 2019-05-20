function filter_unbinned(phase)
import matlab.io.*

% load data directory
files = dir(sprintf('../../data/NH/%s/selected_data/*fit',phase));
file_array=[];

[Nrows, Ncols] = size(files);

% open a txt file to save the name of the unbinned images
filename = fopen(sprintf('../../data/NH/txt/%s_unbinned.txt',phase), 'w');

% count how many unbinned files. Hopefully not too many...
k = 0;

for ifile=1:Nrows
    
    file_to_read= fits.openFile(sprintf('../../data/NH/%s/selected_data/%s',phase,files(ifile).name));
    bin_size_char = fits.readKey(file_to_read, 'NAXIS1');
    bin_size = str2double(bin_size_char);
    
    if bin_size > 256
        fprintf(filename,'%s\n',files(ifile).name);
        k = k + 1;
    end
    fits.closeFile(file_to_read);
    
end

k
return

