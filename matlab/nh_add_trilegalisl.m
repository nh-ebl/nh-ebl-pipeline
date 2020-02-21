function triout=nh_add_trilegalisl(data, paths)
  
  filein = sprintf('%sisl/%s.mat',paths.tridir,data.header.timestamp);

  if numel(dir(filein))
  
    load(filein);

    triout = islout;

  else
    
    disp(sprintf('Trilegal ISL file not available for %s.',data.header.timestamp));
    
    triout.isltotmean = NaN;
    triout.isltoterr = NaN;
    triout.islmaskedmean = NaN;
    triout.islmaskederr = NaN;
    
  end
    
    
end
