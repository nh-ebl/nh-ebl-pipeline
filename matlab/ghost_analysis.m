% This script reads in ghost x/y location, radius, and total position from
% spreadsheet, retrieves x/y location of bright stars just outside FOV from data files, and
% compares to find relationship

% Symons, May 2019

%% Import data from spreadsheet
% Import ghost data from previously made spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /data/symons/nh_data/mat/ghost_info.xlsx
%    Worksheet: 0plan_elon_data_wsol
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2019/05/08 14:36:25

%% Import the data
[~, ~, ghostinfo] = xlsread('/data/symons/nh_data/mat/ghost_info.xlsx','0plan_elon_data_wsol');
ghostinfo = ghostinfo(2:end,:);
ghostinfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),ghostinfo)) = {''};

idx = cellfun(@ischar, ghostinfo);
ghostinfo(idx) = cellfun(@(x) string(x), ghostinfo(idx), 'UniformOutput', false);

%% Clear temporary variables
clearvars idx;

%% Begin ghost analysis
% Import paths for data location
paths = get_paths();

% Load data directory
datafiles = dir(sprintf('%s*.mat',paths.datadir));

% Preallocate space for variables
ghostdist = zeros(size(datafiles,1),1);
ghostpart = zeros(size(datafiles,1),1);
starrad = zeros(size(datafiles,1),1);
starxadj = zeros(size(datafiles,1),1);
staryadj = zeros(size(datafiles,1),1);
starposadj = zeros(size(datafiles,1),1);
ghostposadj = zeros(size(datafiles,1),1);
ghostxadj = zeros(size(datafiles,1),1);
ghostyadj = zeros(size(datafiles,1),1);
starcentdistall = zeros(size(datafiles,1),1);
ghostcentdistall = zeros(size(datafiles,1),1);
starghostdistall = zeros(size(datafiles,1),1);
star2xadj = NaN*ones(size(datafiles,1),1);
star2yadj = NaN*ones(size(datafiles,1),1);
star2posadj = NaN*ones(size(datafiles,1),1);

% For all data files
for ifile=1:size(datafiles,1)
    
    % Print current file number
    disp(sprintf('On file %d of %d.',ifile,size(datafiles,1)));
    
    % Load data
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    % Append ghost x and y position, radius, total position, and if partial
    data.ghost.ghostx = ghostinfo{ifile,9};
    data.ghost.ghosty = ghostinfo{ifile,10};
    data.ghost.ghostrad = ghostinfo{ifile,11};
    data.ghost.ghostdist = ghostinfo{ifile,12};
    data.ghost.ghostpartial = ghostinfo{ifile,13};
    
    % If no ghost, do nothing
    if strcmp(data.ghost.ghostx , '') == 1 || data.ghost.ghostx == 0
        fprintf('No ghost.');
    else
        % If ghost exists, save ghost position
        ghostdist(ifile,1) = data.ghost.ghostdist;
        % If ghost is partial, save that information
        if strcmp(data.ghost.ghostpartial , 'partial ') == 1
            ghostpart(ifile,1) = 1;
            % If ghost not visible in set where other ghosts exist, save that
            % info
        elseif data.ghost.ghostx == 0 && strcmp(data.ghost.ghostpartial , 'partial ') ~= 1
            ghostpart(ifile,1) = 0.5;
        end
        
        % Calculate number of potential bright stars contributing to ghost
        numstars = size(data.ghost.brightmag,2);
        
        % Preallocate space for list of star distance to center, star
        % distance to ghost, and x and y coord of star
        stardistcent = zeros(1,numstars);
        stardistghost = zeros(1,numstars);
        starx = zeros(numstars,1);
        stary = zeros(numstars,1);
        
        % For all stars, retrieve x and y position
        for j = 1:numstars
            
            x = data.ghost.brightxpix(1,j);
            y = data.ghost.brightypix(1,j);
            
            % Adjust coordinates for large (654x654) grid
            if x < 0
                x = 199 + x;
            else
                x = x + 199;
            end
            if y < 0
                y = 199 + y;
            else
                y = y + 199;
            end
            
            % Save star x and y adjusted positions
            starx(j,1) = x;
            stary(j,1) = y;
            
            % Calculate distance from center pixel to star pixel and
            % distance from star pixel to ghost pixel
            stardistcent(1,j) = sqrt((x-(199+128))^2 + (y-(199+128))^2);
            stardistghost(1,j) = sqrt((x-(199+data.ghost.ghostx))^2 + (y-(199+(256-data.ghost.ghosty)))^2);
        end
        
        % Save star to center and star to ghost distances
        data.ghost.stardistcent = stardistcent;
        data.ghost.stardistghost = stardistghost;
        % Find which star (if more than one) is closest to ghost (assuming
        % that is the cause of the ghost)
        [M,I] = min(stardistghost);
        % Save the star/ghost distance and x and y coords for the closest
        % star to ghost
        starrad(ifile,1) = M;
        starxadj(ifile,1) = starx(I,1);
        staryadj(ifile,1) = stary(I,1);
        % For all stars, save star/cent distance, ghost/cent distance, and
        % star/ghost distance
        starcentdistall(ifile,1) = stardistcent(1,I);
        ghostcentdistall(ifile,1) = sqrt(((199+128) - (199+data.ghost.ghostx))^2 + ((199+128) - (199+data.ghost.ghosty))^2);
        starghostdistall(ifile,1) = stardistghost(1,I);
        % Save star total position, ghost total position, and ghost x and y
        % adjusted coords
        starposadj(ifile,1) = sqrt(starxadj(ifile,1)^2 + staryadj(ifile,1)^2);
        ghostposadj(ifile,1) = sqrt((199+data.ghost.ghostx)^2 + (199+data.ghost.ghosty)^2);
        ghostxadj(ifile,1) = 199+data.ghost.ghostx;
        ghostyadj(ifile,1) = 199+data.ghost.ghosty;
        
        % If more than one star per ghost (never seen more than 2), save
        % info for star that's farther away
        if( length(stardistghost) > 1 ) %check out 2 star options instances
            
            % Find star with max distance
            [M,I] = max(stardistghost);
            %             star2rad(ifile,1) = M;
            % Save x and y adjusted coords and adjusted total position
            star2xadj(ifile,1) = starx(I,1);
            star2yadj(ifile,1) = stary(I,1);
            star2posadj(ifile,1) = sqrt(star2xadj(ifile,1)^2 + star2yadj(ifile,1)^2);
            
            % Make plot of ghost location and bright star locations - blue for
            % regular ghost and red for partial ghost
            % this plot is for second star options, just to check them out
            %             h = figure(1);
            %             clf;
            %             x1=200;
            %             x2=455;
            %             y1=200;
            %             y2=455;
            %             xbox = [x1, x2, x2, x1, x1];
            %             ybox = [y1, y1, y2, y2, y1];
            %             plot(xbox, ybox, 'k-', 'LineWidth', 3);
            %             hold on;
            %             xlim([0,654]);
            %             ylim([0,654]);
            %             pbaspect([1 1 1]);
            %             scatter(starx,stary,'y','filled');
            %             scatter(starx(I,1),stary(I,1),'g','filled');
            %             if ghostpart(ifile,1) == 1
            %                 scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
            %             elseif ghostpart(ifile,1) == 0
            %                 scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'b','filled');
            %             end
            %             title(sprintf('%s',data.header.rawfile));
            %             if( ifile < 134 && ifile > 120 )
            %                 pause;
            %             end
        end
        
        % Make plot of ghost location and bright star locations - blue for
        % regular ghost and red for partial ghost
        %         h = figure(1);
        %         set(h,'visible','off');
        %         clf;
        %         x1=200;
        %         x2=455;
        %         y1=200;
        %         y2=455;
        %         xbox = [x1, x2, x2, x1, x1];
        %         ybox = [y1, y1, y2, y2, y1];
        %         plot(xbox, ybox, 'k-', 'LineWidth', 3);
        %         hold on;
        %         xlim([0,654]);
        %         ylim([0,654]);
        %         pbaspect([1 1 1]);
        %         scatter(starx,stary,'y','filled');
        %         if ghostpart(ifile,1) == 1
        %             scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
        %         elseif ghostpart(ifile,1) == 0
        %             scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'b','filled');
        %         end
        %         title(sprintf('%s',data.header.rawfile));
        %         ext = '.png';
        %         imagename = sprintf('%s%s%s',paths.ghostdir,data.header.timestamp,ext);
        %         print(h,imagename, '-dpng');
        
    end
    
    % Save new data to mat files
    %     save(sprintf('%s%s',paths.datadir,datafiles(ifile).name),'data');
    
end

% Plot ghost total position vs distance from center to closest star
figure(2);
scatter(starcentdistall(starcentdistall~=0),ghostposadj(ghostposadj~=0));
xlabel('Distance from center to star');
ylabel('Ghost position');

% Plot distance from center to ghost vs distance from center to star
figure(3);
scatter(starcentdistall(starcentdistall~=0),ghostcentdistall(ghostcentdistall~=0));
xlabel('Distance from center to star');
ylabel('Distance from center to ghost');

% Plot ghost total position vs star total position
figure(4);
scatter(starposadj(starposadj~=0),ghostposadj(ghostposadj~=0));
xlabel('Star position');
ylabel('Ghost position');

% Plot distance from star to ghost vs distance from center to star
figure(5);
scatter(starcentdistall(starcentdistall~=0),starghostdistall(starghostdistall~=0));
xlabel('Distance from center to star');
ylabel('Distance from star to ghost');
% Optionally, add second (farther away) star data to this plot
% hold on;
% scatter(star2posadj,ghostposadj,'MarkerEdgeColor','r');

% Plot star x position vs ghost x position
figure(6);
scatter(starxadj(starxadj~=0), ghostxadj(ghostxadj~=0));
xlabel('Star x position');
ylabel('Ghost x position');

% Plot star y position vs ghost y position
figure(7);
scatter(staryadj(staryadj~=0), ghostyadj(ghostyadj~=0));
xlabel('Star y position');
ylabel('Ghost y position');

drawnow;
