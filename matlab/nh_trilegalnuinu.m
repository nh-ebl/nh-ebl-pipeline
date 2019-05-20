function nh_trilefalnuinu()

  paths = get_paths();
  
  for fieldnum=5:5
    
    load(sprintf('lookup/isltrilegal_%02d.mat',fieldnum));
  
    for ifile=1:10

      mymag = V(ifile).R;
      
      mymag = mymag(mymag > 17.1);
      %mymag = mean([V(ifile).V,V(ifile).R],2);
      
      % convert to flux
      flux = 3e8 .* 1e-26 .* 1e9 .* 3064.*10.^(-mymag/2.5) .* (180./pi).^2 ...
	  ./ 640e-9;
      
      mags = [1:0.5:ceil(max(mymag))+1];
      totflux = zeros(numel(mags)-1,1);
      
      for imag=1:numel(mags)-1
	 whpl = find(mymag < mags(imag));
	 totflux(imag) = sum(flux(whpl));	 
       end

       normflux = totflux ./ max(totflux);
       mags = mags(1:end-1);
       
       if ifile==1
	 figure(1); clf
	 plot(mags,totflux)
	 hold on
	 title(sprintf('Field %d',fieldnum))
	 figure(2); clf
	 plot(mags,normflux)
	 hold on
	 title(sprintf('Field %d',fieldnum))
       else
	 figure(1);
	 plot(mags,totflux);
	 figure(2);
	 plot(mags,normflux);
	 drawnow
       
       end
       
       figure(1);
       pause(1);
       screen2png(sprintf('plots/isl1_%02d.png',fieldnum));

       figure(2);
       pause(1);
       screen2png(sprintf('plots/isl2_%02d.png',fieldnum));

       
     end
    

     
     disp(sprintf('Finished filenum %d',fieldnum));


    
    
    
  end
  
  dbstop
  
  
  
  
% end
