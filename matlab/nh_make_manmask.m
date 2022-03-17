function data= nh_make_manmask(data,paths)

%Plot masked data with optional mouse movement returns value
h = figure(1);
clf;
imagesc(data.data.*(~data.mask.starmask & ~data.mask.linemask & ~data.mask.statmask & ~data.mask.ghostmask))
% imagesc(data.data.*(~data.mask.onemask))
set(h,'visible','off');
% set (gcf, 'WindowButtonMotionFcn', @mouseMove);
a = colorbar;
a.Label.String = 'Intensity [DN]';
pbaspect([1 1 1]);
xlabel('LORRI X Pixels');
ylabel('LORRI Y Pixels');
caxis([-10,10]);
%         title(sprintf('%s',data.header.rawfile));
% grid minor;
% title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
title(sprintf('Field: %d',data.header.fieldnum));
set(gca,'YDir','normal');
ext = '.png';
if not(isfolder([paths.clipmaskdir]))
    mkdir([paths.clipmaskdir])
end
imagename = sprintf('%s%s%s',paths.clipmaskdir,data.header.timestamp,ext);
% imagename = sprintf('%s%s%s',paths.selectmaskdir,data.header.timestamp,ext);
print(h,imagename, '-dpng');

mask = ones(256,256);

data.header.timestamp


if strcmp(data.header.timestamp,'2454001.0080307')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.58.*(jpix+33)+30 & ipix > 0.58.*(jpix+33)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454001.0080770')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.59.*(jpix+43)+30 & ipix > 0.59.*(jpix+43)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454001.0081233')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.62.*(jpix+45)+30 & ipix > 0.62.*(jpix+45)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454001.1468039')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.58.*(jpix+7)+30 & ipix > 0.58.*(jpix+7)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454001.1468502')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.6.*(jpix+9)+30 & ipix > 0.6.*(jpix+9)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454001.1468965')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.6.*(jpix+20)+30 & ipix > 0.6.*(jpix+20)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454110.6986963')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.75.*(jpix+35)+40 & ipix > 0.75.*(jpix+35)-40
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454110.6988236')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.75.*(jpix+45)+40 & ipix > 0.75.*(jpix+45)-40
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454110.6989510')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.75.*(jpix+35)+40 & ipix > 0.75.*(jpix+35)-40
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454110.7299463')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454110.7300736')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454110.7302010')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
                mask(ipix,jpix) = 0;
            end
        end
    end
end
if strcmp(data.header.timestamp,'2454379.5347853')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454379.5349010')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454379.5350168')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454379.5351325')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454379.5352483')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454379.5764520')
    mask(125:165,130:180) = 0;
end
if strcmp(data.header.timestamp,'2454379.5765677')
    mask(125:165,130:180) = 0;
end
if strcmp(data.header.timestamp,'2454379.5766834')
    mask(125:165,130:180) = 0;
end
if strcmp(data.header.timestamp,'2454379.5767992')
    mask(125:165,130:180) = 0;
end
if strcmp(data.header.timestamp,'2454379.5769149')
    mask(125:165,130:180) = 0;
end
if strcmp(data.header.timestamp,'2454379.6181186')
    mask(155:175,130:160) = 0;
end
if strcmp(data.header.timestamp,'2454379.6182344')
    mask(155:175,130:160) = 0;
end
if strcmp(data.header.timestamp,'2454379.6183501')
    mask(155:175,130:160) = 0;
end
if strcmp(data.header.timestamp,'2454379.6184659')
    mask(155:175,130:160) = 0;
end
if strcmp(data.header.timestamp,'2454379.6185816')
    mask(155:175,130:160) = 0;
end
if strcmp(data.header.timestamp,'2454380.1597853')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454380.1599010')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454380.1600168')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454380.1601325')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454380.1602483')
    mask(1:30,1:30) = 0;
end
if strcmp(data.header.timestamp,'2454380.2014520')
    mask(130:160,140:170) = 0;
end
if strcmp(data.header.timestamp,'2454380.2015677')
    mask(130:160,140:170) = 0;
end
if strcmp(data.header.timestamp,'2454380.2016835')
    mask(130:160,140:175) = 0;
end
if strcmp(data.header.timestamp,'2454380.2017992')
    mask(130:160,140:175) = 0;
end
if strcmp(data.header.timestamp,'2454380.2019149')
    mask(130:160,140:175) = 0;
end
if strcmp(data.header.timestamp,'2454755.5105859')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
end
if strcmp(data.header.timestamp,'2454755.5107016')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
end
if strcmp(data.header.timestamp,'2454755.5108174')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
end
if strcmp(data.header.timestamp,'2455370.7917252')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
end
if strcmp(data.header.timestamp,'2455370.7918409')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
end
if strcmp(data.header.timestamp,'2455370.7919567')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
end
if strcmp(data.header.timestamp,'2457583.7727459')
    mask(181:186, 82:95)= 0;
    mask(181:187, 80:83)= 0;
    mask(181:187, 76:81)= 0;
    mask(182:187, 72:77)= 0;
    mask(183:189, 69:73)= 0;
    mask(184:190, 66:70)= 0;
    mask(185:191, 61:67)= 0;
    mask(184:191, 57:60)= 0;
    mask(187:192, 53:58)= 0;
    mask(188:192, 48:54)= 0;
    mask(191:192, 44:47)= 0;
    mask(185:197, 34:39)= 0;
    mask(186:190, 29:30)= 0;
    mask(193:196, 29:34)= 0;
    mask(194:198, 22:28)= 0;
    mask(200:202, 20:22)= 0;
    mask(197:198, 17:20)= 0;
    mask(199:201, 14:17)= 0;
    
end
% if field with quaoar, mask large and small galaxies in field
if strcmp(data.astrom.target,'QUAOAR')
    bigrad = 15; % radius of big galaxy
    smolrad = 10; % radius of small galaxy
    [ypix_biggal, xpix_biggal] = radec2pix(220.7416833,4.8894667, data.astrom); % get pix coords of big galaxy
    [ypix_smolgal, xpix_smolgal] = radec2pix(220.7616667,4.7656667, data.astrom); % get pix coords of small galaxy
    
    % Mask the big galaxy
    %create submask of object (create array of 0's and 1's,
    %where the 1's represnt to location of the objects in the
    %submask). Basically a sqaure of 0's with a circle of 1's
    [rr, cc] = meshgrid(1:2*bigrad+1);
    Circle = sqrt((rr-bigrad-1).^2+(cc-bigrad-1).^2)<=bigrad;
    
    %combined the submask (C) and mask (Z) where xpix, ypix is
    %the center of the object
    for i = 1:(2*bigrad+1)
        xcurr = round(xpix_biggal-bigrad-1+i);
        if xcurr < 1 || xcurr > data.astrom.imagew;
            continue;
        end
        for j = 1:(2*bigrad+1)
            ycurr = round(ypix_biggal-bigrad-1+j);
            if ycurr < 1 || ycurr > data.astrom.imageh;
                continue;
            end
            if Circle(i,j) == 1
                mask(ycurr,xcurr) = 0;
            end
        end
    end
    
    % Mask the small galaxy
    %create submask of object (create array of 0's and 1's,
    %where the 1's represnt to location of the objects in the
    %submask). Basically a sqaure of 0's with a circle of 1's
    [rr, cc] = meshgrid(1:2*smolrad+1);
    Circle = sqrt((rr-smolrad-1).^2+(cc-smolrad-1).^2)<=smolrad;
    
    %combined the submask (C) and mask (Z) where xpix, ypix is
    %the center of the object
    for i = 1:(2*smolrad+1)
        xcurr = round(xpix_smolgal-smolrad-1+i);
        if xcurr < 1 || xcurr > data.astrom.imagew;
            continue;
        end
        for j = 1:(2*smolrad+1)
            ycurr = round(ypix_smolgal-smolrad-1+j);
            if ycurr < 1 || ycurr > data.astrom.imageh;
                continue;
            end
            if Circle(i,j) == 1
                mask(ycurr,xcurr) = 0;
            end
        end
    end
end


% figure(2); clf
% imagesc(mask)
% colorbar
% % caxis([0,2000])

manmask = ~mask;

if not(isfolder([paths.mandir]))
    mkdir([paths.mandir])
end
fileout = sprintf('%s%s.mat',paths.mandir,data.header.timestamp);
save(fileout,'manmask');

end
