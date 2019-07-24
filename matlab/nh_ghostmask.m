%This function takes a data (.mat) file for one LORRI image with an optical ghost
%and creates a mask for the ghost based on known location of bright star in
%range to cause a ghost, previously determined relationships between ghost
%and star position, and previously determined maximum possible ghost radius

%All previously determined values were obtained using the ghost analysis
%script based on all presently available data with optical ghosts

%Symons 2019

function ghostmask = nh_ghostmask(data)

% Calculate number of potential bright stars contributing to ghost
numstars = size(data.ghost.brightmag,2);

% Preallocate space for list of star distance to center and x and y coord of star
stardistcent = zeros(1,numstars);
starx = zeros(numstars,1);
stary = zeros(numstars,1);

% For all stars, retrieve x and y position
for j = 1:numstars
    
    x = data.ghost.brightxpix(1,j);
    y = data.ghost.brightypix(1,j);
    
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
if M <= 268
    
    % Save the star x and y coords for the closest
    % star to ghost and adjust for origin at center of FOV
    starx = data.ghost.brightxpix(1,I) + 199 - 327;
    stary = data.ghost.brightypix(1,I) + 199 - 327;

    % Save the star magnitude
    starmag = data.ghost.brightmag(1,I);

    % Calculate ghost position and magnitude based on star position and
    % magnitude
    % Fit equations were calculated in ghost analysis script
    xghost = 0.1200 * starx + 8.0592;
    yghost = 0.1240 * stary + 1.9207;
    magghost = 1.1473 * starmag + 8.1212;

    % Create mask for ghost
    % Meshgrid based on image size
    [xgrid, ygrid] = meshgrid(1:size(data.data,2), 1:size(data.data,1));
    % Create mask where meshgrid values are less than radius distance from
    % source (adjusting back to origin in lower left corner)
    % Radius is determined by maximum possible ghost radius
    mask = ((xgrid-(xghost+128)).^2 + (ygrid-(yghost+128)).^2) <= 18.^2;
    % Set mask values to 1
    ghostmask = zeros(256);
    ghostmask(mask) = 1;
    
else
    % If star out of range, don't make a ghost mask
    ghostmask = zeros(256);
end

end