%This is a program to retrieve header values from fits files
%It includes the special functionality of being able to take the number of
%targets in the field and loop through and retrieve their names 

%Symons, 2019

clear variables 
close all

%==================User Settings==================
%Linux Path

path = '/data/symons/nh_data_new/';
phase = '';

% path = '/data/symons/NH_old_data/NH/';
% phase = '';


%Windows Path
% path = 'C:\\nh\\NH_old_data\\NH\\regist\\';
% phase = ''; %make it empty to keep it simple

%-----text file-----
%fitsDataFileName = '0FITsDataSummary.txt'; %set the file name for the dataArray of the FITs file keywords (0 makes it at the top!)
%fitsDataFileDelimiter = '\t'; %sets the delimiter between value entries
%-----csv file (Excel opens it nicely! - also csv stands for comma seperated value)-----
fitsDataFileName = '0exp_time_data_new.csv'; %set the file name for the dataArray of the FITs file keywords (0 makes it at the top!)
fitsDataFileDelimiter = ','; %sets the delimiter between value entries

%Make a list of all of the needed keywords in FITs file
fitsKeywordList = {... %start this off
    'SPCBRRA'... %Right Ascension (Boresight RA EME J2000)
    'SPCBRDEC'... %Declination (Boresight DEC EME J2000)
    'SPCUTCAL'... %Date and time in UTC
    'SOL_ELON'... %Solar Elongation
    'EXPTIME'... %Exposure Time
    'TARGET'... %Intended Target
    'TRGFOVN'... %Number of targets in view
    'TRGFOV00'... %MUST BE LAST! Actual target in frame, can be more than one - uses TRGFOV00, TRGFOV01 naming scheme
    }; %end of the keyword list
%Note that keywords not found in the FITs file won't cause errors, just empty places for the data!    


%==================Read the FITs Files==================
%Get the names of the FITs files in the directory
filename = sprintf('%s%s/regist/selected_data/*fit',path,phase); %prep the file name looking system - newest data
% filename = sprintf('%s%s/regist/selected_data/good/*fit',path,phase); %prep the file name looking system - new data
files = dir(filename); %get list of all FITs files in the directory

%Preallocate arrays
dataArray_extraColumns = 3; %extra columns for the data array at the start
dataArray{length(files),length(fitsKeywordList)+dataArray_extraColumns} = []; %preallocate cell array with data (+1 for file name)
%these are used with TRARGFOVN etc.
fitsKeywordListOrig = fitsKeywordList; %record original list
trgNumMax = 1; %start off with a max of 1 targets in FOV
FLG_fitsKeywordListEdited = 0; %flag that denotes that the fits keyword list was edited for more targets

%1st - loop through each FITs file found in the directory
for( i = 1:length(files) )
%     i
    info= fitsinfo(sprintf('%s%s/regist/selected_data/%s',path,phase,files(i).name)); %get FITs file header info - newest data
%     info= fitsinfo(sprintf('%s%s/regist/selected_data/good/%s',path,phase,files(i).name)); %get FITs file header info - new data

    dataArray{i,1} = files(i).name; %record the file name in the first entry
    
    %2nd - loop through each FITs keyword
    j = 1; %while loop used so can change length of fitsKeywordList on fly as needed
    while( j <= length(fitsKeywordList) )
        fitsKeywordIndex =  find(strcmp(info.PrimaryData.Keywords(:,1),fitsKeywordList{j})); %find the index of the keyword
        
        if( ~isempty(fitsKeywordIndex) ) %keeps this from running if the keyword wasn't found
            dataArray(i,j+dataArray_extraColumns) = info.PrimaryData.Keywords(fitsKeywordIndex,2); %record the value of the keyword
            %+dataArray_extraColumns for file name in 1 spot and extra data
        end %end if
        
        %Special check for the number of targets in the FITs file
        if( ~isempty(fitsKeywordIndex) && strcmp(fitsKeywordList{j}, 'TARGFOVN') )
            trgNum = info.PrimaryData.Keywords{fitsKeywordIndex,2}; %get the number out for the current file
            fitsKeywordList{1,length(fitsKeywordListOrig)+trgNum-1} = []; %increase size of the keyword list for the current file
            if( trgNum > trgNumMax )
                dataArray{length(files),size(dataArray,2)+(trgNum - trgNumMax)} = []; %add more cells to hold new cell data
                trgNumMax = trgNum; %update that number to new max
            end
            
            for( k = 2:trgNum ) %increment through the number and add to the keyword list
                fitsKeywordList{j+k-1} = ['TRGFOV',num2str(k-1,'%02.f')]; %create the naming string
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

%3rd - loop through each FITs file found in the directory
loc_ra = find(contains(fitsKeywordList,'SPCBRRA'))+dataArray_extraColumns; %get index of RA for dataArray (+1 for extra data columns)
loc_dec = find(contains(fitsKeywordList,'SPCBRDEC'))+dataArray_extraColumns; %get index of DEC for dataArray (+1 for extra data columns)
for( i = 1:length(files) )
    gacoords = coco([dataArray{i,loc_ra},dataArray{i,loc_dec}],'j2000.0','g','d','d'); %calculate galactic coords
    dataArray{i,2} = gacoords(1);
    dataArray{i,3} = gacoords(2);
end
fitsKeywordList = ['L','B',fitsKeywordList]; %add on new galactic coordinates


%==================Create the Line Formats for the Output File==================
%create the format spec for the header
headerSpec = '%s'; %prep the string
for( j = 1:length(fitsKeywordList) )
    headerSpec = [headerSpec,[fitsDataFileDelimiter,'%s']]; %tack on a %s for each header
end
headerSpec = [headerSpec,'\r\n']; %tack on a line end (\r\n order makes Windows' Notepad happy!)

%create the format spec for the data
formatSpec = '%s'; %prep the string
for( j = 1:length(fitsKeywordList) )
    dataClass = class(dataArray{1,j+1});
    
    if( strcmp(dataClass,'double') == 1 )
        dataAbrev = [fitsDataFileDelimiter,'%f']; % %f means float, which is a double
    elseif( strcmp(dataClass,'char') == 1 )
        dataAbrev = [fitsDataFileDelimiter,'%s']; % %s means string
    else
        error('ERROR: Data class not supported.'); %throw an error that an unsupported data class occured
    end %end if
        
    formatSpec = [formatSpec,dataAbrev]; %tack on the format spec
end
formatSpec = [formatSpec,'\r\n']; %tack on galactic coords & a line end (\r\n order makes Windows' Notepad happy!)


%==================Write the FITs Data Summary File==================
fileID = fopen(fitsDataFileName,'w'); %overwrite current file, if there is one
fprintf(fileID,headerSpec,'fileName',fitsKeywordList{1,:}); %write the header line
for( i = 1:length(files) )
    
    fprintf(fileID,formatSpec,dataArray{i,:}); %write the data lines
    
end %end i


fclose(fileID); %close file




