

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function align_maps.m
%%  This function aligns map2 with the same astrometry and pixelization as 
%%  map1; 
%%
%%  Inputs:
%%  map1 -  The map representing the pixelization and astrometry of the 
%%  resulting map.
%%  map2 -  The map to be aligned.
%%  astrometry1 - The astrometry of the map1.
%%  astrometry2 - The astrometry of the map2.
%%  mask1 - The mask of map1.
%%  mask2 - The mask of map2.
%%
%%  Outputs:
%%  map - The map containing the repixelized map2 that has the 
%%          same astrometry as map1.
%%  mask - the joint mask of map1 and shifted map2.  Pixels where there
%%         is no overlap are masked.
%%  
%%  Changelog:
%%  -v1.0, 19/01/11, Joseph Smidt, original version
%%  -v1.1, Sept 12 2012, MZ, fixed this up to handle masking properly.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [map mask] = align_maps(map1,map2,astrometry1,astrometry2,...
                                 mask1,mask2)

  % Extract astrometry arrays.
  %if numel(map2 ~= 1024.^2)
  
  %    ra1  = astrometry1.ra;
  %    dec1 = astrometry1.dec;
  %    ra2  = astrometry2.ra;
  %    dec2 = astrometry2.dec;
  %else
      
  [X,Y] = meshgrid([1:astrometry1.naxis1],[1:astrometry1.naxis2]);
  
  X = transpose(X);
  Y = transpose(Y);
    
  [ra1,dec1] = pix2radec(astrometry1,X(:),Y(:));
    
  %ra1  = (reshape(ra,astrometry1.naxis1,astrometry1.naxis2));
  %dec1 = (reshape(dec,astrometry1.naxis1,astrometry1.naxis2));

  if astrometry2.naxis1 .* astrometry2.naxis2 ~= 1024.^2
       ra2  = astrometry2.ra;
       dec2 = astrometry2.dec;
  else
      [X,Y] = meshgrid([1:astrometry2.naxis1],[1:astrometry2.naxis2]);
  
      X = transpose(X);
      Y = transpose(Y);
    
      [ra2,dec2] = pix2radec(astrometry2,X(:),Y(:));  
  end
  
  %ra2  = (reshape(ra,astrometry2.naxis1,astrometry2.naxis2));
  %dec2 = (reshape(dec,astrometry2.naxis1,astrometry2.naxis2));
  
  %if size(ra1,1) == 2400 || size(ra2,1) == 2400 
  %dbstop
      %end
  
  % An alternative way to do the same thing. I was just getting griddata
  % being slightly faster on my machine.
  % F = TriScatteredInterp(ra2(:),dec2(:),map2(:),'linear');
  % map = F(ra1,dec1);

  map2p = map2;
  map2p(abs(map2p) > 5e3) = 0;
  
  % Align maps - linear is important since (a) it's the right thing to do
  % and (b) is makes nans in regions where there are no data.
  mapp = griddata(ra2,dec2,map2p,ra1,dec1,'linear');
  
  nmask = ~isfinite(mapp);

  map = griddata(ra2,dec2,map2p,ra1,dec1,'natural');
  map(nmask) = 0;
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % now shift mask2 to match mask 1
  mask2p = round(griddata(ra2,dec2,double(mask2),ra1,dec1,'natural'));
  nmask2p = ~isfinite(mask2p);
  
  mask2p(nmask2p) = 1;
  mask2p = logical(mask2p);

  % and or them together
  mask = mask2p | mask1 | nmask;
    
end
