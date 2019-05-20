function [astrometry] = fits_to_astrometry(astro_head)
%,field,inst,flight)
%funtion from ciber/analysis/Matlab/reuc_astro_align.m 
%orignal author: Mike Zemcov
  
  header = astro_head.PrimaryData.Keywords;
  nentries = length(header);
  
  %astrometry.field = field;
  %astrometry.inst = inst;
  %astrometry.flight = flight;
  astrometry.header = header;
    
  for ihead=1:nentries
    thisentry = lower(header{ihead,1});
    if ~strcmp(thisentry,'comment') & ...
          ~strcmp(thisentry,'history') & ~strcmp(thisentry,'end')
      %Added this if statement, setfield gave error on strings that start with
      %underscore. Since those entries arent used, we can ignore them. 
      if strncmpi(thisentry,'_',1) == 0
      	astrometry = setfield(astrometry,thisentry,header{ihead,2});
      end	
    end
  end
  
  %[X,Y] = meshgrid([1:astrometry.naxis1],[1:astrometry.naxis2]);
  
  %X = transpose(X);
  %Y = transpose(Y);
    
  %[ra,dec] = pix2radec(astrometry,X(:),Y(:));
    
  %astrometry.ra  = (reshape(ra,astrometry.naxis1,astrometry.naxis2));
  %astrometry.dec = (reshape(dec,astrometry.naxis1,astrometry.naxis2));
    
end
