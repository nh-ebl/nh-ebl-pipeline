function nh_darklooker()

  paths = get_paths();

  datafiles = dir(sprintf('%s*.mat',paths.darkdir));

  ndarkfiles = size(datafiles,1);
  
  darkdate = zeros(ndarkfiles,1);
  darkutc = cellstr('');
  darkexp = zeros(ndarkfiles,1);
  darkfield = zeros(ndarkfiles,1);
  darkd = zeros(ndarkfiles,1);
  
  summer = zeros(4,1);
  
  for ifile=1:ndarkfiles
    
    load(sprintf('%s%s',paths.darkdir,datafiles(ifile).name));
    
    darkdate(ifile) = data.header.date_jd-data.header.launch_jd;
    darkutc{ifile} = data.header.date_cal;
    darkexp(ifile) = data.header.exptime;
    darkd(ifile) = data.distance.d_sun_NH;
    
    if darkdate(ifile) < 94
      darkfield(ifile) = 1;
      summer(1) = summer(1) + 1;
    end
    if darkdate(ifile) > 94 & darkdate(ifile) < 96
      darkfield(ifile) = 2;
      summer(2) = summer(2) + 1;
    end
    if darkdate(ifile) > 102 & darkdate(ifile) < 103
      darkfield(ifile) = 3;
      summer(3) = summer(3) + 1;      
    end
    if darkdate(ifile) > 103 & darkdate(ifile) < 104
      darkfield(ifile) = 4;
      summer(4) = summer(4) + 1;
    end
    
    mystring = sprintf('%f, %s, %d, %d, %f',...
	darkdate(ifile),darkutc{ifile},darkfield(ifile),darkexp(ifile),...
	darkd(ifile));
    disp(mystring)
    
  end
  
  dbstop
  