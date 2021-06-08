% Finds variable name used in a folder of .m files. Determines if variable is defined in the .m file,
%   and records the line number the variable was defined at, the file name, 
%   and the number of times the variable was used in the file.
%
% Inputs:
%   varName = variable name in string format to look for
%   folder (Optional) = folder to look for, if not provided the function will use the current working directory
% Outputs:
%   info = structure with 4 fields
%       files = file names of the .m files that use the variable
%       counts = number of times the variable was used each file
%       defined = if the variable was defined in the file, was loaded from a structure, or was just called
%       definedLineNum = line number the variable was defined at

function info = codeTrace(varName,folder)
    if( nargin < 2 )
        sprintf(['Warning: Folder not defined, using current working directory:\n',strrep(pwd,'\','\\')])
        folder = pwd; % Set folder as current working directory if not defined
    end
    
    delimList = {' ','=',',',';','(','[','{',')',']'}; % List of allowed delimiters that signify the end of a variable name
    
    if( ispc )
        files = dir([folder,'\*.m']); %windows gets \
    else
        files = dir([folder,'/*.m']); %linux gets /
    end
    
    files = {files.name}; % Get just the file names
    files = files'; % Flip
    try
        funName = dbstack; % Get this file's name
        funName = funName.file; % Get this file's name
    catch 
        funName = 'codeTrace.m'; % Use hard-coded default if dbstack fails (prob not run as a function)
    end
    funLoc = strfind(files,funName); % Look for this file's name
    funLoc = find(~cellfun('isempty',funLoc)); % Get where this file's name is
    files(funLoc) = []; % Delete the entry - don't search itself
    
    infoCounts = zeros(length(files),1,'int64'); % Preallocate
    [infoDefined{1:length(files)}] = deal('called'); % Prep
    infoDefined = infoDefined'; % Flip
    infoDefinedLineNum = zeros(length(files),1,'int64'); % Preallocate
    
    for i = 1:length(files)
        fID = fopen(files{i});
        text = textscan(fID, '%s', 'Delimiter', '\n'); % Get the text
        fclose(fID);
        text = text{1,1};
%         founds_function = strfind(text,'function'); % Find instances
%         founds_function = find(~cellfun('isempty',founds_function)); % Get non-empty cells
%         logi_function = false; % Prep
%         for j = 1:length(founds_function)
%             if( strfind(text{founds_function(j)},'function') == 1 ) % Find the place
%                 logi_function = true; % Set to true
%             end
%         end
        for deli = 1:length(delimList)
            varNameDelimed = [varName,delimList{deli}]; % Add a delimiter
            founds = strfind(text,varNameDelimed); % Find instances
            foundsLines = find(~cellfun('isempty',founds)); % Get non-empty cells
            for j = 1:length(foundsLines)
                textLine = text{foundsLines(j)}; % Get the text line
                founds_varName = strfind(textLine,varNameDelimed); % Get where the variable is mentioned in the line
                founds_comment = strfind(textLine,'%'); % Get where the comments start in the line
                founds_equals = strfind(textLine,'='); % Get where the equals sign is in the line
                founds_equalsLogi = [strfind(textLine,'=='),strfind(textLine,'>='),strfind(textLine,'<=')]; % Get where the equals sign is in the line
                founds_period = strfind(textLine,'.'); % Get where a period is in the line
                founds_function = strfind(strtrim(textLine),'function '); %Get where function is in the line
                if( ~isempty(founds_function) && (founds_function ~= 1 ) )
                    %function must be the 1st thing for it to be a function defintion
                    founds_function = []; %empty it out, the word is probably in a comment or something
                end
                
                %--- Clear out logical operations that involve = signs ---
                if( ~isempty(founds_equalsLogi) )
                    removr = zeros(length(founds_equals),1,'logical'); % Preallocate
                    for jk = 1:length(founds_equals)
                        removr(jk) = any(abs(founds_equals(jk) - founds_equalsLogi) < 2); % If location differences are less than 2, the equals is part of a logical operation
                    end
                    founds_equals = founds_equals(~removr); % Remove any equals that are logical operators
                end
                
                %--- Remove any entries after comments ---
                if( isempty(founds_comment) == 0 )
                    logi_comment = founds_varName < min(founds_comment); % Only keep variable mentions before the comments begin
                    founds_varName = founds_varName(logi_comment); % Retain only pre-comment stuff
                end
                
                if( ~isempty(founds_varName) ) % Only continue if the entire line wasn't commented out 
                    if( ~isempty(founds_equals) ) % Only check for defining if there's an equals sign involved
                        if( isempty(founds_function) ) % Only continue if the entire line isn't detected as a function definition line
                            %--- Determine if defined in line ---
                            logi_period = any( (founds_varName < founds_equals) ) & ~isempty(founds_period); % Only set defined if it is only used before the equals (and it wasn't loaded from a struct)
                            if( (logi_period == 1) && ~isempty(strfind(infoDefined{i},'called')) )
                                infoDefined{i} = 'loaded from struct'; % Set loaded from struct
                                infoDefinedLineNum(i) = foundsLines(j); % Record the line number
        %                     elseif( (logi_period == 1) && isempty(strfind(infoDefined{i},'called')) )
        %                         infoDefined{i} = [infoDefined{i},' | !!multiple!! loaded from struct']; % Set loaded from struct
        %                         infoDefinedLineNum(i) = foundsLines(j); % Record the line number
                            else
                                if( ~isempty(founds_equals) )
                                    logi_equals = all( (founds_varName < founds_equals) ) & isempty(founds_period); % Only set defined if it is only used before the equals (and it wasn't loaded from a struct)
                                    if( (logi_equals == 1) && ~isempty(strfind(infoDefined{i},'called')) )
                                        infoDefined{i} = 'defined'; % Set defined here
                                        infoDefinedLineNum(i) = foundsLines(j); % Record the line number
        %                             elseif( (logi_equals == 1) && isempty(strfind(infoDefined{i},'called')) )
        %                                 infoDefined{i} = [infoDefined{i},' | !!multiple!! defined']; % Set defined here
        %                                 infoDefinedLineNum(i) = foundsLines(j); % Record the line number
                                    end
                                end
                            end
                        else
                            % If founds_function is not empty, it's a function definiton line
                            founds_varName = strfind(textLine,varName); % Get where the variable is mentioned in the line, no limits
                            if( any(founds_varName < founds_equals) && any(founds_varName > founds_equals) && ~isempty(strfind(infoDefined{i},'called')) )
                                % The variable name is before and after the equals sign, so it's imported, possibly modified, and then exported back out
                                infoDefined{i} = 'function input & output'; % Set function input and output here here
                                infoDefinedLineNum(i) = foundsLines(j); % Record the line number
                            elseif( ~any(founds_varName < founds_equals) && any(founds_varName > founds_equals) && ~isempty(strfind(infoDefined{i},'called')) )
                                % The variable name is after the equals but not before the equals, which means it is a  function input
                                infoDefined{i} = 'function input'; % Set function input here
                                infoDefinedLineNum(i) = foundsLines(j); % Record the line number
                            elseif( any(founds_varName < founds_equals) && ~any(founds_varName > founds_equals) && ~isempty(strfind(infoDefined{i},'called')) )
                                % The variable name is before the equals but not after the equals, which means it s a function output and has to be defined in the function by default
                                infoDefined{i} = 'function output & defined'; % Set function output and defined here
                                infoDefinedLineNum(i) = foundsLines(j); % Record the line number
                            end
                        end
                    end
                    %--- Count occurances ---
                    infoCounts(i) = infoCounts(i) + length(founds_varName); % Count the times the variable name is listed
                end
            end
        end
    end
    
    keepr = infoCounts > 0; % Get what to keep
    files = files(keepr); % Keep only files with counts
    infoCounts = infoCounts(keepr); % Keep only files with counts
    infoDefined = infoDefined(keepr); % Keep only files with counts
    infoDefinedLineNum = infoDefinedLineNum(keepr); % Keep only files with counts
    
    info = struct('files',files,'counts',num2cell(infoCounts),'defined',infoDefined,'definedLineNum',num2cell(infoDefinedLineNum)); % Set the struct to return
end