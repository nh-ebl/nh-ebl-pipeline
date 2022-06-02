%This function retrieves header values from fits files, creates png images
%from those fits files, and appends the retrieved info into the title and
%figure annotation. It includes a custom function to scale the colors of
%the image based on a percentile threshold

function add_info(phase,path)

FLG_overwrite = false; %if set to false will skip already made images, set to true to overwrite

% Rescale axes CLim based on displayed image portion's CData
function rescaleAxesClim(hImage, threshold)
    % Get the displayed image portion's CData
    CData = hImage.CData;
    hAxes = hImage.Parent;
    XLim = fix(hAxes.XLim);
    YLim = fix(hAxes.YLim);
    rows = min(max(min(YLim):max(YLim),1),size(CData,1)); % visible portion
    cols = min(max(min(XLim):max(XLim),1),size(CData,2)); % visible portion
    CData = CData(unique(rows),unique(cols));
    CData = CData(:);  % it's easier to work with a 1d array
 
    % Find the CLims from this displayed portion's CData
    CData = sort(CData(~isnan(CData)));  % or use the Stat Toolbox's prctile()
    thresholdVals = [threshold, 1-threshold];
    thresholdIdxs = fix(numel(CData) .* thresholdVals);
    CLim = CData(thresholdIdxs);
 
    % Update the axes
    hAxes.CLim = CLim;
end

%Make a list of all of the needed keywords in FITs file
fitsKeywordList = {... %start this off
    'EXPTIME'... %Exposure Time
    'SPCUTCAL'... %Spacecraft mid-obs time, UTC [Cal d]
    'TARGET'... %Intended Target
    'TRGFOVN'... %Number of targets in view
    'TRGFOV00'... %MUST BE LAST! Actual target in frame, can be more than one - uses TRGFOV00, TRGFOV01 naming scheme
    }; %end of the keyword list
%Note that keywords not found in the FITs file won't cause errors, just empty places for the data!


%==================Read the FITs Files==================
%Get the names of the FITs files in the directory
fitsfilename = sprintf('%s%s/selected_data/*fit',path,phase); %prep the file name looking system
fitsfiles = dir(fitsfilename); %get list of all FITs files in the directory

%Preallocate arrays
dataArray{length(fitsfiles),length(fitsKeywordList)+1} = []; %preallocate cell array with data (+1 for file name)
%these are used with TRARGFOVN etc.
fitsKeywordListOrig = fitsKeywordList; %record original list
trgNumMax = 1; %start off with a max of 1 targets in FOV
FLG_fitsKeywordListEdited = 0; %flag that denotes that the fits keyword list was edited for more targets

%1st - loop through each FITs file found in the directory
for( i = 1:length(fitsfiles) )
    info= fitsinfo(sprintf('%s%s/selected_data/%s',path,phase,fitsfiles(i).name)); %get FITs file header info
    
    dataArray{i,1} = fitsfiles(i).name; %record the file name in the first entry
    
    %2nd - loop through each FITs keyword
    j = 1; %while loop used so can change length of fitsKeywordList on fly as needed
    while( j <= length(fitsKeywordList) )
        fitsKeywordIndex =  find(strcmp(info.PrimaryData.Keywords(:,1),fitsKeywordList{j})); %find the index of the keyword
        
        if( ~isempty(fitsKeywordIndex) ) %keeps this from running if the keyword wasn't found
            dataArray(i,j+1) = info.PrimaryData.Keywords(fitsKeywordIndex,2); %record the value of the keyword
            %+1 for file name in 1 spot
        end %end if
        
        %Special check for the number of targets in the FITs file
        if( ~isempty(fitsKeywordIndex) && strcmp(fitsKeywordList{j}, 'TRGFOVN') )
            trgNum = info.PrimaryData.Keywords{fitsKeywordIndex,2}; %get the number out for the current file
            fitsKeywordList{1,length(fitsKeywordListOrig)+trgNum-1} = []; %increase size of the keyword list for the current file
            if( trgNum > trgNumMax )
                dataArray{length(fitsfiles),size(dataArray,2)+(trgNum - trgNumMax)} = []; %add more cells to hold new cell data
                trgNumMax = trgNum; %update that number to new max
            end
            
            for( k = 2:trgNum ) %increment through the number and add to the keyword list
                fitsKeywordList{j+k} = ['TRGFOV',num2str(k-1,'%02.f')]; %create the naming string
            end %end for k
            
            FLG_fitsKeywordListEdited = 1; %turn on the edited flag
        end %end if
        
        j = j+1; %increment
    end %end for j
    
    if( FLG_fitsKeywordListEdited == 1 ) %if statement to return the keyword list to its original type
        fitsKeywordList = fitsKeywordListOrig; %overwrite original
        FLG_fitsKeywordListEdited = 0; %reset flag
    end %end if
    
end %end for i

%Iterate through fits files and create pngs
for idate= 1:length(fitsfiles)
    %for saving the image later
    filename = dataArray{idate,1};
    newfilename = strrep(filename,'.fit','.png'); %Replace .fit extension with .png
    imagename = sprintf('%s%s/png_images/%s',path,phase,newfilename);
    
    if( (FLG_overwrite == true) || ~isfile(imagename) )
        fitsfile = fitsfiles(idate).name;%Get file name
        
        data1=fitsread(sprintf('%s%s/selected_data/%s',path,phase,fitsfile)); %Read fits file
        
        %Retrieve exposure time, date of exposure, and target info
        exptime = dataArray{idate, 2};
        date = dataArray{idate, 3};
        trg = dataArray{idate, 4};
        numtrg = dataArray{idate, 5};
        acttrg = '';
        
        %Create string of actual targets in FOV
        for j = 1:numtrg
            acttrg = [acttrg,' ', dataArray{idate,j+5}];
        end
        
        %Create png image
        h = figure(1);
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15.2 11.4]) %Define image resolution
        clf;
        himage = imagesc((data1));
        set(h,'visible','off');
        colormap(flipud(gray));
        colorbar;
        rescaleAxesClim(himage, 0.01)
        %caxis('auto'); %Automatic color scaling
        t2 = annotation('textbox',[0.13,0.075,0,0],'string',acttrg,'FitBoxToText','on'); %Create annotation with target info
        t2.LineStyle = 'none';
        title(['Date Taken: ',date,'; Exp Time: ',num2str(exptime),'; Intended Target: ',trg,'; Num Targets in FOV: ',num2str(numtrg),]); %Put info in title
    
        
        %save as .png file
        if not(isfolder(sprintf('%s%s/png_images/',path,phase)))
            mkdir(sprintf('%s%s/png_images/',path,phase))
        end
        print(h,imagename, '-dpng', '-r100'); %Save as png file with 100x resolution
        
        pause(1);
    end
    
end

return
end


