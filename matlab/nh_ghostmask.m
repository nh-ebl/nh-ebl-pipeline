%This function takes a data (.mat) file for one LORRI image with an optical ghost
%and creates a mask for the ghost based on known location of bright star in
%range to cause a ghost, previously determined relationships between ghost
%and star position, and previously determined maximum possible ghost radius

%All previously determined values were obtained using the ghost analysis
%script based on all presently available data with optical ghosts

%Symons 2019

function [ghostmask,data] = nh_ghostmask(data,paths)

load('run_params.mat','params')
if params.err_mags == 1
    % Choose ghost dir based on err_flag
    ghostdir = data.(params.err_str).ghost;
else
    ghostdir = data.ghost;
end
    
% Calculate number of potential bright stars contributing to ghost
numstars = size(ghostdir.brightmag,1);

% Preallocate space for list of star distance to center and x and y coord of star
stardistcent = zeros(1,numstars);
starx = zeros(numstars,1);
stary = zeros(numstars,1);

% For all stars, retrieve x and y position
for j = 1:numstars
    
    x = ghostdir.brightxpix(j,1);
    y = ghostdir.brightypix(j,1);
    
    % Adjust coordinates for large (654x654) grid
    
    x = 199 + x;
    
    y = 199 + y;
    
    % Calculate distance from center pixel to star pixel and
    % distance from star pixel to ghost pixel
    stardistcent(1,j) = sqrt((x-(199+128))^2 + (y-(199+128))^2);
end

% Find which star (if more than one) is closest to center (assuming
% that is the cause of the ghost)
[M,I] = min(stardistcent);

% If star is in range to cause a ghost, mask ghost
if M <= 268 %18.22 arcmin
    
    % Save the star x and y coords for the closest
    % star to ghost and adjust for origin at center of FOV
    starx = ghostdir.brightxpix(I,1) + 199 - 327;
    stary = ghostdir.brightypix(I,1) + 199 - 327;
    
    % Save the star magnitude
    starmag = ghostdir.brightmag(I,1);
    
    % Calculate ghost position and magnitude based on star position and
    % magnitude
    % Fit equations were calculated and saved in ghost analysis script
    xghost = data.ghost.fitx(1) * starx + data.ghost.fitx(2);
    yghost = data.ghost.fity(1) * stary + data.ghost.fity(2);
    magghost = 1.1473 * starmag + 8.1212;
    
    % Calculate distance from ghost center to center of fov
    ghostdistcent = sqrt((128-(xghost+128))^2 + (128-(yghost+128))^2);
    ghostdir.ghostdistcent = ghostdistcent;
    
    % Create mask for ghost
    % Meshgrid based on image size
    [xgrid, ygrid] = meshgrid(1:size(data.data,2), 1:size(data.data,1));
    % Create mask where meshgrid values are less than radius distance from
    % source (adjusting back to origin in lower left corner)
    % Radius was determined by maximum possible ghost radius (18), now extended
    % to 21.5 to make up for coord prediction having some error
    mask = ((xgrid-(xghost+128)).^2 + (ygrid-(yghost+128)).^2) <= 21.5.^2;
    % Set mask values to 1
    ghostmask = zeros(256);
    ghostmask(mask) = 1;
    
else
    % If star out of range, don't make a ghost mask
    ghostmask = zeros(256);
    ghostdir.ghostdistcent = 0;
end

if params.err_mags == 1
    % Choose ghost dir based on err_flag
    data.(params.err_str).ghost = ghostdir;
else
    data.ghost = ghostdir;
end

end