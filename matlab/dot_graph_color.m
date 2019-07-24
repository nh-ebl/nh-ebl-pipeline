%thayer 06/05/19
%attempt at making dot graph to just show first and last of every group of
%ghosts
clear all
close all

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
starposadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostposadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostxadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostyadj = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starcentdistall = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostcentdistall = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
starghostdistall = zeros((size(datafiles,1)+size(ndatafiles,1)),1);
star2xadj = NaN*ones((size(datafiles,1)+size(ndatafiles,1)),1);
star2yadj = NaN*ones((size(datafiles,1)+size(ndatafiles,1)),1);
star2posadj = NaN*ones((size(datafiles,1)+size(ndatafiles,1)),1);
brightmag= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
ghostmag= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
cts= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
Flux= zeros((size(datafiles,1)+size(ndatafiles,1)),1);
% For old data files
for ifile=1:size(datafiles,1)
    
    % Print current file number
    fprintf('On file %d of %d.\n',ifile,size(datafiles,1));
    
    % Load data
    load(sprintf('%s%s',paths.datadir,datafiles(ifile).name));
    
    % Append ghost x and y position, radius, total position, and if partial
    data.ghost.ghostx = oldghostinfo{ifile,14};
    data.ghost.ghosty = oldghostinfo{ifile,15};
    data.ghost.ghostrad = oldghostinfo{ifile,16};
    data.ghost.ghostdist = oldghostinfo{ifile,17};
    data.ghost.ghostpartial = oldghostinfo{ifile,18};
    [data.ghost.ghostra,data.ghost.ghostdec]=pix2radec(data.astrometry,data.ghost.ghostx,data.ghost.ghosty);
    % If no ghost, do nothing
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
    end
     %
        cts(ifile,1) = ap_photom(data.data.*~data.mask.onemask,data.ghost.ghostx,data.ghost.ghosty,data.ghost.ghostrad,2,3, data, paths);
        Flux= (cts(ifile,1)./data.astrom.exptime);
        ghostmag(ifile,1)=(-2.5*log10(Flux)+20);
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
        brightmag((ifile),1) = data.ghost.brightmag(1,I);
        
        
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
%                         h = figure;
% %                         clf;
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
        % Plot star mag vs ghost mag
  h = figure(2);
      set(h,'visible','on');
  xlim([3.5,6.6]);
  xlabel('Star Magnitude');
  ylim([13,18]);
  ylabel('Ghost Magnitude');
               
                hold on;
                 
                if (brightmag((ifile),1)<4) 
                    gg4= scatter(brightmag((ifile),1),ghostmag((ifile),1),'y','filled');
                    hold on;
                elseif ((brightmag((ifile),1)>4) && (brightmag((ifile),1)<6)) 
                    gg6= scatter(brightmag((ifile),1),ghostmag((ifile),1),'m','filled');
                    hold on;
                elseif ((brightmag((ifile),1)>6) && (brightmag((ifile),1)<6.2)) 
                    gg8= scatter(brightmag((ifile),1),ghostmag((ifile),1),'c','filled');
                    hold on;
                elseif (brightmag((ifile),1)>6.2) 
                    gg10= scatter(brightmag((ifile),1),ghostmag((ifile),1),'g','filled');
                    hold on;
                end
        
  % Make plot of ghost location and bright star locations - blue for
%         % regular ghost and red for partial ghost
                h = figure(1);
                set(h,'visible','on');
%                 clf;
                x1=200;
                x2=455;
                y1=200;
                y2=455;
                xbox = [x1, x2, x2, x1, x1];
                ybox = [y1, y1, y2, y2, y1];
                plot(xbox, ybox, 'k-', 'LineWidth', 3);
                hold on;
                xlim([0,654]);
                ylim([0,654]);
                pbaspect([1 1 1]);
                    
                if ghostpart(ifile,1) == 0 && (data.ghost.ghostx)>=0 && (data.ghost.ghostx)<=116;
                    g1=scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'y','filled');
                end
               if (data.ghost.ghostx)>=0 && (data.ghost.ghostx)<=116;
                   scatter(starx,stary,'y','filled');
               end
               if (data.ghost.ghostx)>=117 && (data.ghost.ghostx)<=160;
                    g2= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'g','filled');
                    scatter(starx,stary,'g','filled');
               end
%                 title(sprintf('%s',data.header.rawfile));
%                 ext = '.png';
%                 imagename = sprintf('%s%s%s',paths.ghostdir,data.header.timestamp,ext);
%                 print(h,imagename, '-dpng');
%              
                
        
       
end

% For new data files
cntr = 7; % counter for figures
for ifile=1:size(ndatafiles,1)
    
    % Print current file number
    disp(sprintf('On file %d of %d.',ifile,size(ndatafiles,1)));
    
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
    end
         %
        cts((ifile+16),1) = ap_photom(data.data.*~data.mask.onemask,data.ghost.ghostx,data.ghost.ghosty,data.ghost.ghostrad,2,3, data, npaths);
        Flux((ifile+16),1)= (cts((ifile+16),1)./data.astrom.exptime);
        ghostmag((ifile+16),1)=(-2.5*log10(Flux((ifile+16),1))+20);
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
        
        %assigning star mag to a variable
        brightmag((ifile+16),1) = data.ghost.brightmag(1,I);
        
        
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
%                         cntr = cntr + 1;
%                         h = figure(cntr);
% %                         clf;
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
        end
%         
  % Plot star mag vs ghost mag
  h = figure(2);
      set(h,'visible','on');
  xlim([3.5,6.6]);
  xlabel('Star Magnitude');
  ylim([13,18]);
  ylabel('Ghost Magnitude');
               
                hold on;
                 
                if (brightmag((ifile+16),1)<4) 
                    gg4= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'y','filled');
                    hold on;
                elseif ((brightmag((ifile+16),1)>4) && (brightmag((ifile+16),1)<6)) 
                    gg6= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'r','filled');
                    hold on;
                elseif ((brightmag((ifile+16),1)>6) && (brightmag((ifile+16),1)<6.2)) 
                    gg8= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'c','filled');
                    hold on;
                elseif (brightmag((ifile+16),1)>6.2) 
                    gg10= scatter(brightmag((ifile+16),1),ghostmag((ifile+16),1),'g','filled');
                    hold on;
                end



%         % Make plot of ghost location and bright star locations - blue for
%         % regular ghost and red for partial ghost
                cntr = cntr + 1;
                h = figure(1);
                set(h,'visible','on');
%                 clf;
                x1=200;
                x2=455;
                y1=200;
                y2=455;
                xbox = [x1, x2, x2, x1, x1];
                ybox = [y1, y1, y2, y2, y1];
                plot(xbox, ybox, 'k-', 'LineWidth', 3);
                xlim([0,654]);
                ylim([0,654]);
                pbaspect([1 1 1]);
                hold on;
                 if ghostpart((ifile+16),1) == 1
%                     g1= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'y','filled');
%                     hold on;
%                     scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'y','filled');
                 elseif (brightmag((ifile+16),1)<6)
                    g3= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'r','filled');
                   
                    scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'r','filled');
                 elseif (brightmag((ifile+16),1)>6)
                    g4= scatter(199+data.ghost.ghostx,199+(256-data.ghost.ghosty),'c','filled');
                 
                    scatter(starxadj((ifile+16),1),staryadj((ifile+16),1),'c','filled');
%             
                 end    
                 
%                 title(sprintf('%s',data.header.rawfile));
%                 ext = '.png';
%                 imagename = sprintf('%s%s%s',npaths.ghostdir,data.header.timestamp,ext);
%                 print(h,imagename, '-dpng');
        
end
    

fit= polyfit(brightmag(brightmag~=0),ghostmag(ghostmag~=0),1);
starfit=linspace(min(brightmag(brightmag~=0)),max(brightmag(brightmag~=0)));
ghostfit=(fit(1)*starfit + fit(2));
plot (starfit, ghostfit);
text(4.5,17,'y=0.6806x+11.6902');
legend([gg4,gg6,gg8,gg10],{'Mag 3.95','Mag 5.79','Mag 6.15','Mag 6.35',});

    % Save new data to mat files
%          save(sprintf('%s%s',npaths.datadir,ndatafiles(ifile).name),'data');
% 



legend([g1,g2,g3,g4,],{'PrePluto Star Mag 3.9555','PrePluto Mag 6.3535','Mag < 5.7967','Mag 6.1574'});

drawnow;




