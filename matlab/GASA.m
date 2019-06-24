% original: ghost_analysis.m
% This script reads in ghost x/y location, radius, and total position from
% excel spreadsheet, retrieves x/y location of bright stars just outside FOV from data files,
% and compares to find relationships
% Symons, May 2019
% These relationships involve x/y locations of ghosts and stars
% all graphs are color coded so they all line up too
% Edits: Thayer June 2019
clear all
close all

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
%% Import the data old
[~, ~, oldghostinfo] = xlsread('/data/symons/NH_old_data/mat/good/old_ghost_info_test.xlsx','old_ghost_info');
oldghostinfo = oldghostinfo(2:end,:);
oldghostinfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),oldghostinfo)) = {''};
idx = cellfun(@ischar, oldghostinfo);
oldghostinfo(idx) = cellfun(@(x) string(x), oldghostinfo(idx), 'UniformOutput', false);

%% Import the data new
[~, ~, ghostinfo] = xlsread('/data/symons/nh_data/mat/ghost_info_test.xlsx','0plan_elon_data_wsol');
ghostinfo = ghostinfo(2:end,:);
ghostinfo(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),ghostinfo)) = {''};
idx = cellfun(@ischar, ghostinfo);
ghostinfo(idx) = cellfun(@(x) string(x), ghostinfo(idx), 'UniformOutput', false);

%% Clear temporary variables
clearvars idx;

%% Begin ghost analysis
% Import paths for data location
paths = get_paths();
npaths=get_paths_new();

% Load both data directories
%old
datafiles = dir(sprintf('%s*.mat',paths.datadir));
%new
ndatafiles= dir(sprintf('%s*.mat',npaths.datadir));

% Preallocate space for variables
ghostdist = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostpart = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starrad = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starxadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
staryadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starxcent=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starycent=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starposadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostposadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostxadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostyadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostxcent=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostycent=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starcentdistall = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostcentdistall = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starghostdistall = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
star2xadj = NaN*ones((size(datafiles,1)+size(ndatafiles,1)),1);
star2yadj = NaN*ones((size(datafiles,1)+size(ndatafiles,1)),1);
star2posadj = NaN*ones((size(datafiles,1)+size(ndatafiles,1)),1);
brightmag= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostmag= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
cts= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostxguess=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostyguess=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostmagguess=zeros((size(datafiles,1)+size(ndatafiles,1)),1);
Flux= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostrad= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
% For old data files
for ifile=1:size(datafiles,1)
    
    % Print current file number
    fprintf('On file %d of %d.\n',ifile,size(datafiles,1));
    
    % Load data
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    % Append ghost x and y position, radius, total position, and if partial
    %these come from the excel files where the locations were predetermined
    data.ghost.ghostx = oldghostinfo{ifile,14};
    data.ghost.ghosty = oldghostinfo{ifile,15};
    data.ghost.ghostrad = oldghostinfo{ifile,16};
    data.ghost.ghostdist = oldghostinfo{ifile,17};
    data.ghost.ghostpartial = oldghostinfo{ifile,18};
    [data.ghost.ghostra,data.ghost.ghostdec]=pix2radec(data.astrometry,data.ghost.ghostx,data.ghost.ghosty);
    
    % If no ghost (nothing in x/y location), do nothing
    if  strcmp(data.ghost.ghostx, 0) == 1
        fprintf('No ghost ');
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
        
        % ghost brightness analysis
        %using ap_photom, the brightness of the ghost is determined and
        %assigned to variables
        cts(ifile,1) = ap_photom(data.data.*~data.mask.onemask,data.ghost.ghostx,data.ghost.ghosty,data.ghost.ghostrad,2,3,data,paths);
        Flux(ifile,1)= (cts(ifile,1)./data.astrom.exptime);
        ghostmag(ifile,1)=(-2.5*log10(Flux(ifile,1))+20);
       
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
        
        %assigning star mag to a variable
        brightmag(ifile,1) = data.ghost.brightmag(1,I);
        ghostrad(ifile,1)= data.ghost.ghostrad;
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
        % testing the location of ghosts using our prediction equation
        ghostxguess(ifile,1)= 0.1208*starxadj(ifile,1)+296.7564;
        ghostyguess(ifile,1)=0.1127*staryadj(ifile,1)+293.9804;
        ghostmagguess(ifile,1)=0.6930*brightmag(ifile,1)+11.5544;
        
        ghostxcent(ifile,1)= ghostxadj(ifile,1)-327;
        ghostycent(ifile,1)= ghostyadj(ifile,1)-327;
        starxcent(ifile,1)= starxadj(ifile,1)-327;
        starycent(ifile,1)= staryadj(ifile,1)-327;
       
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
%                         h = figure;
%                         clf;
%                         x1=200;
%                         x2=455;
%                         y1=200;
%                         y2=455;
%                         xbox = [x1, x2, x2, x1, x1];
%                         ybox = [y1, y1, y2, y2, y1];
%                         plot(xbox, ybox, 'k-', 'LineWidth', 3);
%                         hold on;
%                         xlim([0,654]);
%                         ylim([0,654]);
%                         pbaspect([1 1 1]);
%                         scatter(starx,stary,'y','filled');
%                         scatter(starx(I,1),stary(I,1),'g','filled');
%                         if ghostpart(ifile,1) == 1
%                             scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
%                         elseif ghostpart(ifile,1) == 0
%                             scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'b','filled');
%                         end
%                         title(sprintf('%s',data.header.rawfile));
%                         if( ifile < 134 && ifile > 120 )
%                             pause;
%                         end
        end
%         
%  Plot star magnitude vs ghost magnitude: Vertical plot; this displays the
%  magnitudes organized by color, with the color representing different
%  star magnitudes
%                 h = figure(3);
%                 set(h,'visible','on');
%                 xlim([3.5,6.6]);
%                 xlabel('Star Magnitude');
%                 ylim([12,16.5]);
%                 ylabel('Ghost Magnitude');
% 
%                 if (brightmag((ifile),1)<4) 
%                     g1= scatter(brightmag(ifile,1),ghostmag(ifile,1),'r','filled');
%                     hold on;
%                 elseif ((brightmag((ifile),1)>4) && (brightmag((ifile),1)<6)) 
%                     g2= scatter(brightmag(ifile,1), ghostmag(ifile,1),'m','filled');
%                     hold on;
%                 elseif ((brightmag((ifile),1)>6) && (brightmag((ifile),1)<6.2)) 
%                     g3= scatter(brightmag(ifile,1), ghostmag(ifile,1),'c','filled');
%                     hold on;
%                 elseif (brightmag((ifile),1)>6.2) 
%                     g4= scatter(brightmag(ifile,1),ghostmag(ifile,1),'b','filled');
%                     hold on;
%                 end
% Dot graph! Inside the box is the ghost location, outside is star location
% but is color coded by the star magnitude
%                 h = figure(2);
%                 set(h,'visible','on');
%                 x1=200;
%                 x2=455;
%                 y1=200;
%                 y2=455;
%                 xbox = [x1, x2, x2, x1, x1];
%                 ybox = [y1, y1, y2, y2, y1];
%                 plot(xbox, ybox, 'k-', 'LineWidth', 3);
%                 xlim([0,654]);
%                 ylim([0,654]);
%                 pbaspect([1 1 1]);
%                 hold on;
%                  if ghostpart((ifile),1) == 1
% %                     gg1= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'y','filled');
% %                     hold on;
% %                     scatter(starxadj(ifile,1),staryadj(ifile,1),'y','filled');
% %                  elseif (brightmag((ifile),1)~=0)
%                      
%                  elseif (brightmag((ifile),1)<4)
%                     g1= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
%                    hold on;
%                     scatter(starxadj((ifile),1),staryadj((ifile),1),'r','filled');
%                     
%                  elseif ((brightmag((ifile),1)>4) && (brightmag((ifile),1)<6))
%                     g2= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'m','filled');
%                  hold on;
%                     scatter(starxadj((ifile),1),staryadj((ifile),1),'m','filled');
%                     
%                  elseif ((brightmag((ifile),1)>6) && (brightmag((ifile),1)<6.2))
%                      g3= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'c','filled');
%                  hold on;
%                     scatter(starxadj((ifile),1),staryadj((ifile),1),'c','filled');
%                     
%                  elseif ((brightmag((ifile),1)>6.2))
%                      g4= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'b','filled');
%                  hold on;
%                     scatter(starxadj((ifile),1),staryadj((ifile),1),'b','filled');
%                  end
                 
%                 title(sprintf('%s',data.header.rawfile));
%                 ext = '.png';
%                 imagename = sprintf('%s%s%s',npaths.ghostdir,data.header.timestamp,ext);
%                 print(h,imagename, '-dpng');
% Save new data to mat files
%          save(sprintf('%s%s',npaths.datadir,ndatafiles(ifile).name),'data');
% %   
if (brightmag(ifile,1)<4)
    gm1= ghostmag(ifile,1);
    gr1= ghostrad(ifile,1);
elseif ((brightmag((ifile),1)>4) && (brightmag((ifile),1)<6))
    gm2= ghostmag(ifile,1);
    gr2= ghostrad(ifile,1);
elseif ((brightmag((ifile),1)>6) && (brightmag((ifile),1)<6.2))
    gm3= ghostmag(ifile,1);
    gr3= ghostrad(ifile,1);
elseif ((brightmag((ifile),1)>6.2))
    gm4= ghostmag(ifile,1);
    gr4= ghostrad(ifile,1);
end

        figure(4);
        hold on;
        xlabel('Ghost Radius');
        xlim([11.5 20]);
        ylabel('Ghost Magnitude');
        ylim([13.5 17]);
        if (brightmag((ifile),1)<4)
            g1= scatter(gr1(gr1~=0),gm1(gm1~=0),'r', 'filled');
        elseif ((brightmag((ifile),1)>4) && (brightmag((ifile),1)<6))
            g2= scatter(gr2,gm2,'m', 'filled');
        elseif ((brightmag((ifile),1)>6) && (brightmag((ifile),1)<6.2))
            g3= scatter(gr3,gm3,'c', 'filled');
        elseif ((brightmag((ifile),1)>6.2))
            g4= scatter(gr4,gm4,'b', 'filled');
        end
        

    end
    


end

% For new data files

for ifile=1:size(ndatafiles,1)
    
    % Print current file number
    fprintf('On file %d of %d.\n',ifile,size(ndatafiles,1));
    
    % Load data
    load(sprintf('%s%s',npaths.datadir,ndatafiles(ifile).name));
    
    % Append ghost x and y position, radius, total position, and if partial
    data.ghost.ghostx = ghostinfo{ifile,14};
    data.ghost.ghosty = ghostinfo{ifile,15};
    data.ghost.ghostrad = ghostinfo{ifile,16};
    data.ghost.ghostdist = ghostinfo{ifile,17};
    data.ghost.ghostpartial = ghostinfo{ifile,18};
    
    % If no ghost, do nothing
    if strcmp(data.ghost.ghostx , '') == 1 || data.ghost.ghostx == 0
        
        fprintf('No ghost.');
    else
        [data.ghost.ghostra,data.ghost.ghostdec]=pix2radec(data.astrometry,data.ghost.ghostx,data.ghost.ghosty);
        % If ghost exists, save ghost position
        ghostdist((ifile+16),1) = data.ghost.ghostdist;
        % If ghost is partial, save that information
        if strcmp(data.ghost.ghostpartial , 'partial ') == 1
            ghostpart((ifile+16),1) = 1;
            % If ghost not visible in set where other ghosts exist, save that
            % info
            
        elseif data.ghost.ghostx == 0 && strcmp(data.ghost.ghostpartial , 'partial ') ~= 1
            ghostpart((ifile+16),1) = 0.5;
       
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
        starrad((ifile+16),1) = M;
        starxadj((ifile+16),1) = starx(I,1);
        staryadj((ifile+16),1) = stary(I,1);

        % For all stars, save star/cent distance, ghost/cent distance, and
        % star/ghost distance
        starcentdistall((ifile+16),1) = stardistcent(1,I);
        ghostcentdistall((ifile+16),1) = sqrt(((199+128) - (199+data.ghost.ghostx))^2 + ((199+128) - (199+data.ghost.ghosty))^2);
        starghostdistall((ifile+16),1) = stardistghost(1,I);
        % Save star total position, ghost total position, and ghost x and y
        % adjusted coords
        starposadj((ifile+16),1) = sqrt(starxadj((ifile+16),1)^2 + staryadj((ifile+16),1)^2);
        ghostposadj((ifile+16),1) = sqrt((199+data.ghost.ghostx)^2 + (199+data.ghost.ghosty)^2);
        ghostxadj((ifile+16),1) = 199+data.ghost.ghostx;
        ghostyadj((ifile+16),1) = 199+data.ghost.ghosty;
        
        ghostxguess((ifile+16),1)= 0.1208*starxadj((ifile+16),1)+296.7564;
        ghostyguess((ifile+16),1)=0.1127*staryadj((ifile+16),1)+293.9804;
        %assigning star magnitude to a variable
        if strcmp(data.ghost.ghostpartial , 'partial ') == 1
            
        elseif data.ghost.ghostx == 0 && strcmp(data.ghost.ghostpartial , 'partial ') ~= 1
            
        else 
        %using ap_photom, the brightness of the ghost is determined and
        %assigned to variables
        cts((ifile+16),1) = ap_photom(data.data.*~data.mask.onemask,data.ghost.ghostx,data.ghost.ghosty,data.ghost.ghostrad,2,3,data,npaths);
        Flux((ifile+16),1)= (cts((ifile+16),1)./data.astrom.exptime);
        ghostmag((ifile+16),1)=(-2.5*log10(Flux((ifile+16),1))+20);
        brightmag((ifile+16),1) = data.ghost.brightmag(1,I);
        ghostmagguess((ifile+16),1)=0.6930*brightmag((ifile+16),1)+11.5544;
        ghostrad((ifile+16),1)= data.ghost.ghostrad;
        
        end
        ghostxcent((ifile+16),1)= ghostxadj((ifile+16),1)-327;
        ghostycent((ifile+16),1)= ghostyadj((ifile+16),1)-327;
        starxcent((ifile+16),1)= starxadj((ifile+16),1)-327;
        starycent((ifile+16),1)= staryadj((ifile+16),1)-327;
        

        % If more than one star per ghost (never seen more than 2), save
        % info for star that's farther away
        if( length(stardistghost) > 1 ) %check out 2 star options instances
            
            % Find star with max distance
            [M,I] = max(stardistghost);
            %             star2rad(ifile,1) = M;
            % Save x and y adjusted coords and adjusted total position
            star2xadj((ifile+16),1) = starx(I,1);
            star2yadj((ifile+16),1) = stary(I,1);
            star2posadj((ifile+16),1) = sqrt(star2xadj((ifile+16),1)^2 + star2yadj((ifile+16),1)^2);
            
            
   
            % Make plot of ghost location and bright star locations - blue for
            % regular ghost and red for partial ghost
            % this plot is for second star options, just to check them out
%                         cntr= cntr+1;
%                         h = figure(cntr);
%                         clf;
%                         x1=200;
%                         x2=455;
%                         y1=200;
%                         y2=455;
%                         xbox = [x1, x2, x2, x1, x1];
%                         ybox = [y1, y1, y2, y2, y1];
%                         plot(xbox, ybox, 'k-', 'LineWidth', 3);
%                         hold on;
%                         xlim([0,654]);
%                         ylim([0,654]);
%                         pbaspect([1 1 1]);
%                         scatter(starx,stary,'y','filled');
%                         scatter(starx(I,1),stary(I,1),'g','filled');
%                         if ghostpart((ifile+16),1) == 1
%                             scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
%                         elseif ghostpart((ifile+16),1) == 0
%                             scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'b','filled');
%                         end
%                         title(sprintf('%s',data.header.rawfile));
%                         if( (ifile+16) < 134 && (ifile+16) > 120 )
%                             pause;
%                         end
%hold off;
        end
%         
% Plot star magnitude vs ghost magnitude: Vertical plot; this displays the
% magnitudes organized by color, with the color representing different
% star magnitudes  
%                h = figure(3);
%                set(h,'visible','on');
%                xlim([3.5,6.6]);
%                xlabel('Star Magnitude');
%                ylim([12,16.5]);
%                ylabel('Ghost Magnitude');
%                
%                 if (brightmag((ifile+16),1)<4) 
%                     g1= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'r','filled');
%                     hold on;
%                 elseif ((brightmag((ifile+16),1)>4) && (brightmag((ifile+16),1)<6)) 
%                     g2= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'m','filled');
%                     hold on;
%                 elseif ((brightmag((ifile+16),1)>6) && (brightmag((ifile+16),1)<6.2)) 
%                     g3= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'c','filled');
%                     hold on;
%                 elseif (brightmag((ifile+16),1)>6.2) 
%                     g4= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'b','filled');
%                     hold on;
%                 end

% Dot graph! Inside the box is the ghost location, outside is star location
% but is color coded by the star magnitude
%                 h = figure(2);
%                 set(h,'visible','on');
%                 x1=200;
%                 x2=455;
%                 y1=200;
%                 y2=455;
%                 xbox = [x1, x2, x2, x1, x1];
%                 ybox = [y1, y1, y2, y2, y1];
%                 plot(xbox, ybox, 'k-', 'LineWidth', 3);
%                 xlim([0,654]);
%                 ylim([0,654]);
%                 pbaspect([1 1 1]);
%                 hold on;
%                  if ghostpart((ifile+16),1) == 1
% %                     gg1= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'y','filled');
% %                     hold on;
% %                     scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'y','filled');
%                  elseif (brightmag((ifile+16),1)<4)
%                     g1= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
%                    hold on;
%                     scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'r','filled');
%                  elseif ((brightmag((ifile+16),1)>4) && (brightmag((ifile+16),1)<6))
%                     g2= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'m','filled');
%                  hold on;
%                     scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'m','filled');
%                  elseif ((brightmag((ifile+16),1)>6) && (brightmag((ifile+16),1)<6.2))
%                     g3= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'c','filled');
%                  hold on;
%                     scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'c','filled');
%                  elseif (brightmag((ifile+16),1)>6.2)
%                     g4= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'b','filled');
%                  hold on;
%                     scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'b','filled');
%                  end
%                  
%                 title(sprintf('%s',data.header.rawfile));
%                 ext = '.png';
%                 imagename = sprintf('%s%s%s',npaths.ghostdir,data.header.timestamp,ext);
%                 print(h,imagename, '-dpng');
% Save new data to mat files
%          save(sprintf('%s%s',npaths.datadir,ndatafiles(ifile).name),'data');

if (brightmag((ifile+16),1)<4)
    gm1= ghostmag((ifile+16),1);
    gr1= ghostrad((ifile+16),1);
elseif ((brightmag((ifile+16),1)>4) && (brightmag((ifile+16),1)<6))
    gm2= ghostmag((ifile+16),1);
    gr2= ghostrad((ifile+16),1);
elseif ((brightmag((ifile+16),1)>6) && (brightmag((ifile+16),1)<6.2))
    gm3= ghostmag((ifile+16),1);
    gr3= ghostrad((ifile+16),1);
elseif ((brightmag((ifile+16),1)>6.2))
    gm4= ghostmag((ifile+16),1);
    gr4= ghostrad((ifile+16),1);
end

        figure(4);
        hold on;
        xlabel('Ghost Radius');
        xlim([11.5 20]);
        ylabel('Ghost Magnitude');
        ylim([13.5 17]);
        if (brightmag((ifile),1)<4)
            g1= scatter(gr1(gr1~=0),gm1(gm1~=0),'r', 'filled');
        elseif ((brightmag((ifile),1)>4) && (brightmag((ifile),1)<6))
            g2= scatter(gr2,gm2,'m', 'filled');
        elseif ((brightmag((ifile),1)>6) && (brightmag((ifile),1)<6.2))
            g3= scatter(gr3,gm3,'c', 'filled');
        elseif ((brightmag((ifile),1)>6.2))
            g4= scatter(gr4,gm4,'b', 'filled');
        end
      
    end
    
end

% These are for the plots inside the for loop.
% The fit is for the magnitude plot(figure 3).
% fit= polyfit(brightmag(brightmag~=0),ghostmag(ghostmag~=0),1);
% starfit=linspace(min(brightmag(brightmag~=0)),max(brightmag(brightmag~=0)));
% ghostfit=(fit(1)*starfit + fit(2));
% figure(3);
% hold on;
% plot (starfit, ghostfit);
% text(4.5,13,'y=0.9931x+8.8923');
legend([g1,g2,g3,g4],{'Mag 3.9555','Mag 5.7967','Mag 6.1574','Mag 6.3535',});
% hold off;
% figure(2);
% hold on;
% legend([g1,g2,g3,g4],{'Mag 3.9555','Mag 5.7967','Mag 6.1574','Mag 6.3535',});
% hold off;

%% Graphs
cntr = 5; % counter for figures
% Plot ghost total position vs distance from center to closest star
% figure(cntr);
% scatter(starcentdistall(starcentdistall~=0),ghostposadj(ghostposadj~=0));
% xlabel('Distance from center to star');
% ylabel('Ghost position');
% 
% % Plot distance from center to ghost vs distance from center to star
% figure(cntr+1);
% scatter(starcentdistall(starcentdistall~=0),ghostcentdistall(ghostcentdistall~=0));
% xlabel('Distance from center to star');
% ylabel('Distance from center to ghost');
% 
% % Plot ghost total position vs star total position
% figure(cntr+2);
% scatter(starposadj(starposadj~=0),ghostposadj(ghostposadj~=0));
% % scatter(starposadjpol,ghostposadjpol);
% xlabel('Star position');
% ylabel('Ghost position');
% 
% % Plot distance from star to ghost vs distance from center to star
% figure(cntr+3);
% scatter(starcentdistall(starcentdistall~=0),starghostdistall(starghostdistall~=0));
% xlabel('Distance from center to star');
% ylabel('Distance from star to ghost');
% % Optionally, add second (farther away) star data to this plot
% % hold on;
% % scatter(star2posadj,ghostposadj,'MarkerEdgeColor','r');

%% Possible polar graphs below ()
%(theta,rho)
[startadj,starradj]=cart2pol(starxcent, starycent);
startadj= rad2deg(startadj);
[ghosttadj,ghostradj]=cart2pol(ghostxcent, ghostycent);
ghosttadj= rad2deg(ghosttadj);

% Plot star r position vs ghost r position
% figure(cntr+4);
% hold on;
% scatter(starradj(starradj~=0), ghostradj(ghostradj~=0));
% fitr=polyfit(starradj(starradj~=0),ghostradj(ghostradj~=0),1);
% starrfit= linspace(min(starradj(starradj~=0)),max(starradj(starradj~=0)));
% ghostrfit= (fitr(1).*starrfit + fitr(2));
% plot(starrfit,ghostrfit);
% xlabel('Star r position');
% ylabel('Ghost r position');
% % polarscatter(starradj(starradj~=0), ghostradj(ghostradj~=0));
% title('[starradj, ghostradj]');
% 
% % Plot star t position vs ghost t position
% figure(cntr+5);
% hold on;
% scatter(startadj(startadj~=0), ghosttadj(ghosttadj~=0));
% xlabel('Star theta position');
% ylabel('Ghost theta position');
% fitt=polyfit(startadj(startadj~=0),ghosttadj(ghosttadj~=0),1);
% startfit= linspace(min(startadj(startadj~=0)),max(startadj(startadj~=0)));
% ghosttfit= (fitt(1).*startfit + fitt(2));
% plot(startfit,ghosttfit);
% text(0,100,'y=0.8886x-5.4691');
% % polarscatter(startadj(startadj~=0), ghosttadj(ghosttadj~=0));
% title('[startadj,ghosttadj]');

% Plot star mag vs ghost mag
% figure(cntr+6);
% scatter(brightmag(brightmag~=0), ghostmag(ghostmag~=0));
% xlabel('Star Magnitude');
% ylabel('Ghost Magnitude');
% xlim([3.9,6.6]);
% ylim([12,16.5]);
% hold on;
% fitmag= polyfit(brightmag(brightmag~=0),ghostmag(ghostmag~=0),1);
% starfit=linspace(min(brightmag(brightmag~=0)),max(brightmag(brightmag~=0)));
% ghostfit=(fitmag(1)*starfit + fitmag(2));
% plot (starfit, ghostfit);
% text(4.5,13,'y=0.9931x+8.8923')
% hold off;

%% Ghost Test Location
%%%re-edit

%  %ghost location x test
%  figure(cntr+7);
%  hold on;
% %  scatter(ghostxadj(ghostxadj~=0),ghostxguess(ghostxguess~=0));
%  fitgx=polyfit(ghostxadj(ghostxadj~=0),ghostxguess(ghostxguess~=0),1);
%  gxeq= (fitgx(1)*ghostxadj(ghostxadj~=0) +fitgx(2));
%  guessmeqx= (ghostxguess(ghostxguess~=0)- gxeq);
% %  scatter(ghostxadj(ghostxadj~=0),guessmeqx);
%      % to plot in polar coordinates use the following three lines
% [ghostxguesspol,gxeqpol]= cart2pol(ghostxguess(ghostxguess~=0), gxeq); 
% [ghostxadjpol,guessmeqxpol]= cart2pol(ghostxadj(ghostxadj~=0),guessmeqx);
% fitgx=polyfit(ghostxadj(ghostxadj~=0),ghostxguess(ghostxguess~=0),1);
% gxeq= (fitgx(1)*ghostxadj(ghostxadj~=0) +fitgx(2));
% guessmeqx= (ghostxguess(ghostxguess~=0)- gxeq);
% scatter(ghostxadjpol,guessmeqxpol);
%  xlabel('xadj');
%  ylabel('guess-fit');
%  hold off;
%  
%  %ghost location y test
%  figure(cntr+8);
%  hold on;
% %  scatter(ghostyadj(ghostyadj~=0),ghostyguess(ghostyguess~=0));
%  fitgy=polyfit(ghostyadj(ghostyadj~=0),ghostyguess(ghostyguess~=0),1);
%  gyeq= (fitgy(1)*ghostyadj(ghostyadj~=0) +fitgy(2));
%  guessmeqy= (ghostyguess(ghostyguess~=0) - gyeq);
% %  scatter(ghostyadj(ghostyadj~=0),guessmeqy,'b');
%     % to plot in polar coordinates use the following two lines
% [ghostyadjpol,guessmeqypol]=cart2pol(ghostyadj(ghostyadj~=0),guessmeqy);
% scatter(ghostyadjpol,guessmeqypol);
%  xlabel('yadj');
%  ylabel('guess-fit');
%  hold off;
%  
%  %ghost location mag test
%  figure(cntr+9);
%  hold on;
% %  scatter(ghostmag,ghostmagguess);
%  fitgmag=polyfit(ghostmag,ghostmagguess,1);
%  gmageq= (fitgmag(1)*ghostmag +fitgmag(2));
%  guessmeq= (ghostmagguess - gmageq);
% %  scatter(ghostmag,guessmeq);
%     % to plot in polar coordinates use the following two lines
% [ghostmagpol,guessmeqpol]=cart2pol(ghostmag,guessmeq);
% scatter(ghostmagpol,guessmeqpol);
%  xlabel('Ghost Mag');
%  ylabel('Guess-fit');
%  hold off;

%% 

%plot ghost magnitude against ghost radius

