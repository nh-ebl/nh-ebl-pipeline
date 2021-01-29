% """
% This function regrids an array by a decimal shift through mathematical attempts
% inputs:
%     inArray - input array of any size
%     shiftXY - (shiftX, shiftY) decimal shift - can be any size, any decimal, no fixegers required
%     padarrayVal *UNUSED* - default 0, can set the value it padarrays with on the edges to a constant number
%     FLG_origSize - default 0, 0 means that the outArray will NOT have the original inArray size (as the array spreads out when shifting)
%                     1 means that the outArray will have the original inArray size using RegridderZen to readjust back (RegridderZen is flux conserving as well)
% @author - RD && TS 2019
% """

function outArray = RegridderShift(inArray, shiftXY, padarrayVal, FLG_origSize)
    inSize = size(inArray); %get the size of the in array
    shiftX = shiftXY(1); %get the X shift
    shiftXpix = fix(ceil(abs(shiftX))); %number of pixels to shift in X by
    if( shiftXpix == 0 )
        shiftXdir = fix(0); %0 for zero
    else
        shiftXdir = fix(round(shiftX/abs(shiftX))); %1 for positive, -1 for negative
        if( shiftXdir == -1 ) %fix so I don't have to write negative shift code directly
            shiftX = rem(shiftX,1); %get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
            if( shiftX < 0 )
                shiftX = 1 + shiftX; %fix it up so it's the positive mirror
            end
        else
            shiftX = rem(shiftX,1); %get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
        end
    end
    shiftY = shiftXY(2); %get the Y shift
    %fix for wrong direction Y shift (was shifting down for +values)
    %shiftY = -shiftY %fixes it
    shiftYpix = fix(ceil(abs(shiftY))); %number of pixels to shift in X by
    if( shiftYpix == 0 )
        shiftYdir = fix(0); %0 for zero
    else
        shiftYdir = fix(round(shiftY/abs(shiftY))); %1 for positive, -1 for negative
        if( shiftYdir == -1 ) %fix so I don't have to write negative shift code directly
            shiftY = rem(shiftY,1); %get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
            if( shiftY < 0 )
                shiftY = 1 + shiftY; %fix it up so it's the positive mirror
            end
        else
            shiftY = rem(shiftY,1); %get just the decimal portion (if a shift of like 2.2 occured, it'll do 2.2-2=0.2)
        end
    end
    
    outSize = [ fix(inSize(1) + abs(shiftYdir)) , fix(inSize(2) + abs(shiftXdir)) ]; %took a bit to make nopython happy
    %tuples are annoying
    
    %create the matrix to shift the X
    if( shiftX ~= 0 )
        shiftMatX = diag( repmat((1-shiftX) , [outSize(2),1]), 0 ); %diagonal of (1-shiftX)
        shiftMatX = shiftMatX + diag( repmat(shiftX , [inSize(2),1]), 1  );  %offset diagonal of shiftX
        shiftMatX = shiftMatX(1:end-1,:);  %cut off the bottom row
        %basically this creates a matrix that combines two values in X - just how the shift does
    else
        shiftMatX = diag( ones([inSize(2),1]), 0  ); %I matrix so nothing changes, since shiftX == 0
    end
    
    %create the matrix to shift the Y
    if( shiftY ~= 0 )
        shiftMatY = diag( repmat(shiftY , [outSize(1),1]), 0  ); %diagonal of shiftY
        shiftMatY = shiftMatY + diag( repmat((1-shiftY) , [inSize(1),1]), -1  ); %offset diagonal of (1-shiftY)
        shiftMatY = shiftMatY(:,1:end-1);  %cut off the right column
        %basically this creates a matrix that combines two values in Y - just how the shift does
    else
        shiftMatY = diag( ones([inSize(1),1]), 0  ); %I matrix so nothing changes, since shiftY == 0
    end
    
    %NOTE THAT padarrayVAL IS UNUSED - ASSUMED TO BE 0. Will need to figure out how to implement if not 0.
    outArray = shiftMatY*inArray*shiftMatX; %do the shift calc in one shot
    
    %now after all of that is donezo, we're gonna shift by whole pixels
    if( ((shiftX ~= 0) && ((shiftXpix-1) > 0)) && (shiftXdir == 1) ) %there's a whole pixel shift in the positive direction (->)
        outArray = padarray(outArray, [ 0 , 0 ; (shiftXpix-1) , 0 ] , padarrayVal); %padarray with 0's on the edge that's needed (left gets padarrayded with 0's since the whole thing moved right (+) by X pixels)
    elseif( ((shiftX ~= 0) && ((shiftXpix-1) > 0)) && (shiftXdir == -1) ) %there's a whole pixel shift in the negative direction (<-)
        outArray = padarray(outArray, [ 0 , 0 ; 0 , (shiftXpix-1) ] , padarrayVal); %padarray with 0's on the edge that's needed (right gets padarrayded with 0's since the whole thing moved left (-) by X pixels)
    end
    if( ((shiftY ~= 0) && ((shiftYpix-1) > 0)) && (shiftYdir == 1) ) %there's a whole pixel shift in the positive direction (->)
        outArray = padarray(outArray, [ 0 , (shiftYpix-1) ; 0 , 0 ] , padarrayVal); %padarray with 0's on the edge that's needed (down gets padarrayded with 0's since the whole thing moved up (+) by X pixels)
    elseif( ((shiftY ~= 0) && ((shiftYpix-1) > 0)) && (shiftYdir == -1) ) %there's a whole pixel shift in the negative direction (<-)
        outArray = padarray(outArray, [ (shiftYpix-1) , 0 ; 0 , 0 ] , padarrayVal); %padarray with 0's on the edge that's needed (up gets padarrayded with 0's since the whole thing moved down (-) by X pixels)
    end
    %special case for 1 pixel shift, maybe I could do it better but this is what I got :3
    if( ((shiftX == 0) && (shiftXpix > 0)) && (shiftXdir == 1) ) %there's a whole pixel shift in the positive direction (->)
        outArray = padarray(outArray, [ 0 , 0 ; (shiftXpix) , 0 ] , padarrayVal); %padarray with 0's on the edge that's needed (left gets padarrayded with 0's since the whole thing moved right (+) by X pixels)
    elseif( ((shiftX == 0) && (shiftXpix > 0)) && (shiftXdir == -1) ) %there's a whole pixel shift in the negative direction (<-)
        outArray = padarray(outArray, [ 0 , 0 ; 0 , (shiftXpix) ] , padarrayVal); %padarray with 0's on the edge that's needed (right gets padarrayded with 0's since the whole thing moved left (-) by X pixels)
    end
    if( ((shiftY == 0) && (shiftYpix > 0)) && (shiftYdir == 1) ) %there's a whole pixel shift in the positive direction (->)
        outArray = padarray(outArray, [ 0 , (shiftYpix) ; 0 , 0 ] , padarrayVal); %padarray with 0's on the edge that's needed (down gets padarrayded with 0's since the whole thing moved up (+) by X pixels)
    elseif( ((shiftY == 0) && (shiftYpix > 0)) && (shiftYdir == -1) ) %there's a whole pixel shift in the negative direction (<-)
        outArray = padarray(outArray, [ (shiftYpix) , 0 ; 0 , 0 ] , padarrayVal); %padarray with 0's on the edge that's needed (up gets padarrayded with 0's since the whole thing moved down (-) by X pixels)
    end
    
    if( FLG_origSize == 1 ) %if original size flag is on, keep the original size fixact
        outArray = RegridderZen(outArray, inSize ); %regrid for the original input size
    end

%display("SUM CHECK!!\nSUM IN: "+num2str(sum(inArray))+"\nSUM OUT: "+num2str(sum(outArray)))
    
end