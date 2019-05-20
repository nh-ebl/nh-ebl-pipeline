% MATLAB script to construct star mask in NH data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C.Nguyen
% 11/22/15 Ver 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
% data = path to data directory
% ndate = number of files in data directory
% file_array = array of the file names
% 
% x_array = x-axis for the re-constructed image
% y_array = y-axis for the re-constructed image
% Z = array of MASKED image
% z = array used to store original image
% Z_image = catalyst to pass image data into
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mask(phase,path)
% load matrix with date, file names, and median and deviation
load(sprintf('%smatrix/%s_med_dev.mat',path,phase));
import matlab.io.*


files = file_array;
med = median_array;
sigma = stdev_array;

% define x- and y-axis for the images. Each image is 256x256 binned pixels
x_array=linspace(0,256,256);
y_array=linspace(0,256,256);

% define z& Z arrays. Z array will store observed data, except for the
% masked region where Z(i,ii) = 0. z array will store only observed data
Z=[];
z=[];

[Nrows,Ncols] = size(files);

%for idate= 1:numel(files) % loop over all selected files
 for idate = 1:Ncols
    % pass the .fit name to name
    name = sprintf('%s%s/selected_data/%s',path,phase,files(idate).name);    
    %destination = sprintf('../../data/NH/%s/masked',phase);
    %copyfile(name,destination);
    
    masked_data = sprintf('%s%s/masked/masked_%s',path,phase, files(idate).name);

    
    %fits.openFile(masked_data)
    % pass image data into Z_image array
    Z_image = single(fitsread(name));
    
    z_cutoff_plus = med(idate) + 2.0*sigma(idate);
    z_cutoff_minus = med(idate) - 2.0*sigma(idate);
    
   
    
    % construct map
    for i=1:256
    
      ii = 1;
    
      for ii=1:256

          
          % re-construct image
          z(i,ii) = Z_image(i,ii);
          %Z(i,ii) = z(i,ii);
          
          % filter out the center we don't want
          if z(i,ii) < z_cutoff_plus && z(i,ii) > z_cutoff_minus
            Z(i,ii) = Z_image(i,ii);
          else
            Z(i,ii) = NaN;
          end
        
         ii = ii+1;
      end
    
      i = i+1;
    end
    
    imagedata = single(Z);
    fitswrite(imagedata, masked_data);
    
    
end

    
    figure(1)
    clf; % clear fig window, so that the window will close after each image is saved
    h = imagesc(x_array,y_array,Z);
    set(h,'alphadata',~isnan(Z));
    colormap(flipud(gray));
    caxis([-5,20]);
    colorbar; % display the colorbar on the images
    
%     figure(2)
%     imagesc(x_array,y_array,Z);
%     colormap(flipud(gray));
%     caxis([0,500]);
%     colorbar; 


return


