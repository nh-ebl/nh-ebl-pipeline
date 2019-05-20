function nh_make_manmask(filenum)

  paths = get_paths();
   
  %datafiles = dir(sprintf('%s*.mat',paths.datadir));

  load(sprintf('%s%s.mat',paths.datadir,filenum));
  
  figure(1); clf
  imagesc(data.image.calimage.*~data.mask.onemask)
  colorbar
    
  mask = ones(256,256);
  
  
  if strcmp(filenum,'2454001.0080307')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.58.*(jpix+33)+30 & ipix > 0.58.*(jpix+33)-30
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end
  if strcmp(filenum,'2454001.0080770')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.59.*(jpix+43)+30 & ipix > 0.59.*(jpix+43)-30
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end
  if strcmp(filenum,'2454001.0081233')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.62.*(jpix+45)+30 & ipix > 0.62.*(jpix+45)-30
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end
  if strcmp(filenum,'2454001.1468039')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.58.*(jpix+7)+30 & ipix > 0.58.*(jpix+7)-30
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end
  if strcmp(filenum,'2454001.1468502')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.6.*(jpix+9)+30 & ipix > 0.6.*(jpix+9)-30
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end
  if strcmp(filenum,'2454001.1468965')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.6.*(jpix+20)+30 & ipix > 0.6.*(jpix+20)-30
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end  
  if strcmp(filenum,'2454110.6986963')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.75.*(jpix+35)+40 & ipix > 0.75.*(jpix+35)-40
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end  
  if strcmp(filenum,'2454110.6988236')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.75.*(jpix+45)+40 & ipix > 0.75.*(jpix+45)-40
	  mask(ipix,jpix) = 0;
	end
      end
    end  
  end
  if strcmp(filenum,'2454110.6989510')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.75.*(jpix+35)+40 & ipix > 0.75.*(jpix+35)-40
	  mask(ipix,jpix) = 0;
	end
      end
    end    
  end
  if strcmp(filenum,'2454110.7299463')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
	  mask(ipix,jpix) = 0;
	end
      end
    end 
  end
  if strcmp(filenum,'2454110.7300736')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
	  mask(ipix,jpix) = 0;
	end
      end
    end 
  end
  if strcmp(filenum,'2454110.7302010')
    for ipix=1:256
      for jpix=1:256
	if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
	  mask(ipix,jpix) = 0;
	end
      end
    end 
  end
  if strcmp(filenum,'2454379.5347853')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454379.5349010')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454379.5350168')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454379.5351325')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454379.5352483')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454379.5764520')
    mask(125:165,130:180) = 0;
  end
  if strcmp(filenum,'2454379.5765677')
    mask(125:165,130:180) = 0;
  end
  if strcmp(filenum,'2454379.5766834')
    mask(125:165,130:180) = 0;
  end
  if strcmp(filenum,'2454379.5767992')
    mask(125:165,130:180) = 0;
  end
  if strcmp(filenum,'2454379.5769149')
    mask(125:165,130:180) = 0;
  end
  if strcmp(filenum,'2454379.6181186')
    mask(155:175,130:160) = 0;
  end
  if strcmp(filenum,'2454379.6182344')
    mask(155:175,130:160) = 0;
  end
  if strcmp(filenum,'2454379.6183501')
    mask(155:175,130:160) = 0;
  end
  if strcmp(filenum,'2454379.6184659')
    mask(155:175,130:160) = 0;
  end
  if strcmp(filenum,'2454379.6185816')
    mask(155:175,130:160) = 0;
  end
  if strcmp(filenum,'2454380.1597853')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454380.1599010')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454380.1600168')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454380.1601325')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454380.1602483')
    mask(1:30,1:30) = 0;
  end
  if strcmp(filenum,'2454380.2014520')
    mask(130:160,140:170) = 0;
  end
  if strcmp(filenum,'2454380.2015677')
    mask(130:160,140:170) = 0;
  end
  if strcmp(filenum,'2454380.2016835')
    mask(130:160,140:175) = 0;
  end
  if strcmp(filenum,'2454380.2017992')
    mask(130:160,140:175) = 0;
  end
  if strcmp(filenum,'2454380.2019149')
    mask(130:160,140:175) = 0;
  end
  if strcmp(filenum,'2454755.5105859')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
  end
  if strcmp(filenum,'2454755.5107016')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
  end
  if strcmp(filenum,'2454755.5108174')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
  end
  if strcmp(filenum,'2455370.7917252')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
  end
  if strcmp(filenum,'2455370.7918409')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
  end
  if strcmp(filenum,'2455370.7919567')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
  end
  
  
  
  figure(2); clf
  imagesc(data.image.calimage.*~data.mask.onemask.*mask)
  colorbar
  caxis([0,2000])
  
  manmask = ~mask;
  
  fileout = sprintf('%s%s.mat',paths.mandir,filenum);
  save(fileout,'manmask');
  
% end
