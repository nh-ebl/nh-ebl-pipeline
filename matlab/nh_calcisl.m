function data=nh_calcisl(data, paths)

  usnoisl = nh_usnoisl(data, paths);
   
  data.isl.usnotot = usnoisl.isltot;
  data.isl.usnototim = usnoisl.totimage;
  data.isl.usnowing = usnoisl.islwing;
  data.isl.usnowingim = usnoisl.wingimage;
  data.isl.usnofaint = usnoisl.islfaint;    
  
%   gsciiisl = nh_gsciiisl(data);
%   
%   data.isl.gsciitot = gsciiisl.isltot;
%   data.isl.gsciitotim = gsciiisl.totimage;
%   data.isl.gsciiwing = gsciiisl.islwing;
%   data.isl.gsciiwingim = gsciiisl.wingimage;
%   data.isl.gsciifaint = gsciiisl.islfaint;   
      
  triout = nh_add_trilegalisl(data, paths);
  
  data.isl.tritotmean = triout.isltotmean;
  data.isl.tritoterr = triout.isltoterr;
  data.isl.trimean = triout.islmaskedmean;
  data.isl.trierr = triout.islmaskederr;

 end