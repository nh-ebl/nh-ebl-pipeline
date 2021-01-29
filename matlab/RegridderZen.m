%
%This function regrids an array to a new array size through mathematical perfection
%inputs:
%    inArray - input array of any size
%    outSize - desired outsize of the array
%@author - RD & TS 2019
%
function outArray = RegridderZen(inArray, outSize)
    inSize = size(inArray); %get the size of the in array
    xFact = inSize(1)/outSize(1); %get the y resize factor
    yFact = inSize(2)/outSize(2); %get the x resize factor
    %outArray = zeros(outSize); %create output array
    
    %========== Y MATRIX CREATION - YES IT USES ALL X STUFF (maybe I named it wrong at the start lol?) ==========
    xRange = 0:xFact:inSize(1)+xFact; %get a range of pixel slices to make
    xRangeAdj = xRange; %copy cause python memory stuff
    for i = 1:length(xRangeAdj) %run through all of the x ranges and adjust as needed
        if( xRangeAdj(i) == fix(xRangeAdj(i)) ) %catch true integers
            xRangeAdj(i) = xRangeAdj(i) - .1; %force any int values to int64(floor to the previous integer for counting reasons
        end
    end
    
    if( inSize(1) ~= outSize(1) )
        yMatrix = zeros( [outSize(1),inSize(1)] );
        for i = 1:1:outSize(1)
            xPixelsInvolved = fix(xRangeAdj(i)):1:fix(xRangeAdj(i+1)); %get the pixels involved
            xPixelPercent = ones(size(xPixelsInvolved)); %get the percentages involved
            if( length(xPixelsInvolved) == 1 ) %if size is one, then need the difference between the xRange(2) and xRange(1) to cover the divide
                xPixelPercent(1) = xRange(i+1) - xRange(i); %sets the only value to be the difference between the percentages
            else
                if( i == 1 ) %if i == 1 it means that xRange(1) = 0, so special +1 added
                    xPixelPercent(1) = ceil(xRange(i)) + 1 - xRange(i); %sets the first value to be the first percentage
                else
                    xPixelPercent(1) = ceil(xRange(i)) - xRange(i); %sets the first value to be the first percentage
                end
                xPixelPercent(end) = xRange(i+1) - fix(xRangeAdj(i+1)); %sets the last value to be the last percentage
                %middle numbers are 1 since full pixels are taken, and if there's only 1 pixel involved the last percentage (one that matters )
            end
            for j = 1:length(xPixelsInvolved)
                yMatrix(i,xPixelsInvolved(j)+1) = xPixelPercent(j); %put in the pixel percent (the extra magic is using xPixelsInvolved as an index)
            end %FOR j
        end %FOR i
    else
        yMatrix = diag( ones([inSize(1),1]), 0  ); %I matrix so nothing changes, since shiftX == 0
    end
    
    
    %========== X MATRIX CREATION - YES IT USES ALL Y STUFF (maybe I named it wrong at the start lol?) ==========
    yRange = 0:yFact:inSize(2); %get a range of pixel slices to make
    yRangeAdj = yRange; %copy cause python memory stuff
    for i = 1:length(yRangeAdj) %run through all of the y ranges and adjust as needed
        if( yRangeAdj(i) == fix(yRangeAdj(i)) ) %catch true integers
            yRangeAdj(i) = yRangeAdj(i) - .1; %force any int values to int64(floor to the previous integer for counting reasons
        end
    end
    
    if( inSize(2) ~= outSize(2) )
        xMatrix = zeros( [inSize(2),outSize(2)] ); %preallocate
        for i = 1:outSize(2)
            yPixelsInvolved = fix(yRangeAdj(i)):1:fix(yRangeAdj(i+1)); %get the pixels involved
            yPixelPercent = ones(size(yPixelsInvolved)); %get the percentages involved
            if( length(yPixelsInvolved) == 1 ) %if size is one, then need the difference between the yRange(2) and yRange(1) to cover the divide
                yPixelPercent(1) = yRange(i+1) - yRange(i); %sets the only value to be the difference between the percentages
            else
                if( i == 1 ) %if j == 1 it means that yRange(1) = 1, so special +1 added
                    yPixelPercent(1) = ceil(yRange(i)) + 1 - yRange(i); %sets the first value to be the first percentage
                else
                    yPixelPercent(1) = ceil(yRange(i)) - yRange(i); %sets the first value to be the first percentage
                end
                 yPixelPercent(end) = yRange(i+1) - fix(yRangeAdj(i+1)); %sets the last value to be the last percentage
                %middle numbers are 1 since full pixels are taken, and if there's only 1 pixel involved the last percentage (one that matters )
            end
            for j = 1:length(yPixelsInvolved)
                xMatrix(yPixelsInvolved(j)+1,i) = yPixelPercent(j); %put in the pixel percent (the extra magic is using yPixelsInvolved as an index)
            end %FOR j
        end %FOR i
        
    else
        xMatrix = diag( ones([inSize(2),1]), 0  ); %I matrix so nothing changes, since shiftX == 0
    end   
    
    %Convert to sparse matricies for speed (mostly 0's)
    if( (inSize(1)>1000) || (inSize(2)>1000) || (outSize(1)>1000) || (outSize(2)>1000) )
        %if array large, use sparse arrays - big speed up (but gotta be big to be worth it)
        yMatrix = sparse(yMatrix); %the resize will use sparse matrix math to use these right, 3.5x speed up on large array
        xMatrix = sparse(xMatrix);
        %========== THE RESIZE ==========
        outArray = yMatrix*inArray*xMatrix; %yMatrix and xMatrix are sparse matricies - do sparse math to go fast
    else
        %smol array and regular, don't use sparse
        %========== THE RESIZE ==========
        % if (all(xMatrix == yMatrix) == False):  % prevents the math if a (size1,size1) -> (size1,size1) was requested [so nothing changes and xMatrix/yMatrix = Identity Matrix]
        outArray = yMatrix*inArray*xMatrix;  % one shot calc the resize
        % else:
        %     outArray = inArray %would be same, just set it to be same
        % end [commented out assuming time saving didn't matter much if rarely will have samesize->samesize]
        % outArray = yMatrix@inArray@xMatrix %one shot calc the resize [matmul may be faster]
    end
end