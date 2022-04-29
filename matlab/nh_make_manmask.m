function data= nh_make_manmask(data,paths)

manmask_mode = 'inter'; %'man' for manually typing in coordinates, 'inter' for interactive coordinate selection
manmask_mode_interOrientation = 'horiz'; %horiz | vert - vert makes the 2 plots appear vertically, horiz makes them appear horizontally

%Plot masked data with optional mouse movement returns value
% h = figure(1);
% clf;
% imagesc(data.data.*(~data.mask.starmask & ~data.mask.linemask & ~data.mask.statmask & ~data.mask.ghostmask))
% % imagesc(data.data.*(~data.mask.onemask))
% set(h,'visible','off');
% % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
% a = colorbar;
% a.Label.String = 'Intensity [DN]';
% pbaspect([1 1 1]);
% xlabel('LORRI X Pixels');
% ylabel('LORRI Y Pixels');
% caxis([-10,10]);
% %         title(sprintf('%s',data.header.rawfile));
% % grid minor;
% % title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
% title(sprintf('Field: %d',data.header.fieldnum));
% set(gca,'YDir','normal');
% ext = '.png';
% if not(isfolder([paths.clipmaskdir]))
%     mkdir([paths.clipmaskdir])
% end
% imagename = sprintf('%s%s%s',paths.clipmaskdir,data.header.timestamp,ext);
% % imagename = sprintf('%s%s%s',paths.selectmaskdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

if( strcmp(manmask_mode,'man') )
    %Plot side-by-side masked data without and with clip mask with optional mouse movement returns value
    h = figure(1);
    clf;
    set(h,'Position',[1150 200 768 1024]); %set fig size [leftOffset bottomOffset horizontalFigSize verticalFigSize]
    ax1 = subplot(2,1,1);
    imagesc(data.data.*(~data.mask.onemask))
    set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    % set(h,'visible','off');
    a = colorbar;
    a.Label.String = 'Intensity [DN]';
    pbaspect([1 1 1]);
    xlabel('LORRI X Pixels');
    ylabel('LORRI Y Pixels');
    caxis([-10,10]);
    set(gca,'YDir','normal');
    th1 = title(''); %place holder
    titlePos = get(th1, 'position');
    titlePosNew = titlePos + [-60, 0, 0];
    set(th1,'Position',titlePosNew);
    set(th1, 'units', 'normal');
    
    ax2 = subplot(2,1,2);
    imagesc(data.data.*(~data.mask.starmask & ~data.mask.linemask & ~data.mask.statmask & ~data.mask.ghostmask & ~data.mask.manmask))
    set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    a = colorbar;
    a.Label.String = 'Intensity [DN]';
    pbaspect([1 1 1]);
    xlabel('LORRI X Pixels');
    ylabel('LORRI Y Pixels');
    caxis([-10,10]);
    %         title(sprintf('%s',data.header.rawfile));
    % grid minor;
    % title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
    % title(sprintf('Field: %d',data.header.timestamp));
    set(gca,'YDir','normal');
    th2 = title(''); %place holder
    titlePos = get(th2, 'position');
    titlePosNew = titlePos + [-60, 0, 0];
    set(th2,'Position',titlePosNew);
    set(th2, 'units', 'normal');

    %link axes for syncronized zoom
    linkaxes([ax1, ax2],'xy')
    %disp(''); %break point catch

    % ext = '.png';
    % if not(isfolder([paths.clipmaskdir]))
    %     mkdir([paths.clipmaskdir])
    % end
    % imagename = sprintf('%s%s%s',paths.clipmaskdir,data.header.timestamp,ext);
    % imagename = sprintf('%s%s%s',paths.selectmaskdir,data.header.timestamp,ext);
    % print(h,imagename, '-dpng');
end

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
elseif strcmp(data.header.timestamp,'2454001.0080770')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.59.*(jpix+43)+30 & ipix > 0.59.*(jpix+43)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454001.0081233')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.62.*(jpix+45)+30 & ipix > 0.62.*(jpix+45)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454001.1468039')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.58.*(jpix+7)+30 & ipix > 0.58.*(jpix+7)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454001.1468502')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.6.*(jpix+9)+30 & ipix > 0.6.*(jpix+9)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454001.1468965')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.6.*(jpix+20)+30 & ipix > 0.6.*(jpix+20)-30
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454110.6986963')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.75.*(jpix+35)+40 & ipix > 0.75.*(jpix+35)-40
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454110.6988236')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.75.*(jpix+45)+40 & ipix > 0.75.*(jpix+45)-40
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454110.6989510')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.75.*(jpix+35)+40 & ipix > 0.75.*(jpix+35)-40
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454110.7299463')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454110.7300736')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454110.7302010')
    for ipix=1:256
        for jpix=1:256
            if ipix < 0.64.*(jpix+65) + 40 & ipix > 0.64.*(jpix+65) - 40
                mask(ipix,jpix) = 0;
            end
        end
    end
elseif strcmp(data.header.timestamp,'2454379.5347853')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454379.5349010')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454379.5350168')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454379.5351325')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454379.5352483')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454379.5764520')
    mask(125:165,130:180) = 0;
elseif strcmp(data.header.timestamp,'2454379.5765677')
    mask(125:165,130:180) = 0;
elseif strcmp(data.header.timestamp,'2454379.5766834')
    mask(125:165,130:180) = 0;
elseif strcmp(data.header.timestamp,'2454379.5767992')
    mask(125:165,130:180) = 0;
elseif strcmp(data.header.timestamp,'2454379.5769149')
    mask(125:165,130:180) = 0;
elseif strcmp(data.header.timestamp,'2454379.6181186')
    mask(155:175,130:160) = 0;
elseif strcmp(data.header.timestamp,'2454379.6182344')
    mask(155:175,130:160) = 0;
elseif strcmp(data.header.timestamp,'2454379.6183501')
    mask(155:175,130:160) = 0;
elseif strcmp(data.header.timestamp,'2454379.6184659')
    mask(155:175,130:160) = 0;
elseif strcmp(data.header.timestamp,'2454379.6185816')
    mask(155:175,130:160) = 0;
elseif strcmp(data.header.timestamp,'2454380.1597853')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454380.1599010')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454380.1600168')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454380.1601325')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454380.1602483')
    mask(1:30,1:30) = 0;
elseif strcmp(data.header.timestamp,'2454380.2014520')
    mask(130:160,140:170) = 0;
elseif strcmp(data.header.timestamp,'2454380.2015677')
    mask(130:160,140:170) = 0;
elseif strcmp(data.header.timestamp,'2454380.2016835')
    mask(130:160,140:175) = 0;
elseif strcmp(data.header.timestamp,'2454380.2017992')
    mask(130:160,140:175) = 0;
elseif strcmp(data.header.timestamp,'2454380.2019149')
    mask(130:160,140:175) = 0;
elseif strcmp(data.header.timestamp,'2454755.5105859')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
elseif strcmp(data.header.timestamp,'2454755.5107016')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
elseif strcmp(data.header.timestamp,'2454755.5108174')
    mask(85:135,100:150) = 0;
    mask(:,119:123) = 0;
    mask(112:114,:) = 0;
elseif strcmp(data.header.timestamp,'2455370.7917252')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
elseif strcmp(data.header.timestamp,'2455370.7918409')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
elseif strcmp(data.header.timestamp,'2455370.7919567')
    mask(109:159,111:161) = 0;
    mask(:,135:137) = 0;
    mask(133:135,:) = 0;
elseif strcmp(data.header.timestamp,'2457583.7727459')
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

% !!!Start of Lauer fields!!!
elseif strcmp(data.header.timestamp,'2458351.2085071')
    mask(156:165, 70:79)= 0;
    mask(164:169, 96:101)= 0;
    mask(32, 55:65)= 0;
    mask(52, 187:193)= 0;
    mask(144:146, 155:162)= 0;
    mask(164:166, 166:171)= 0;
    mask(101:104, 155:160)= 0;
    mask(164:166, 166:171)= 0;
    mask(32, 24)= 0;
    mask(34, 39)= 0;
    mask(36, 56)= 0;
    mask(19:20, 65)= 0;
    mask(20:21, 35:37)= 0;
    mask(39:41, 222:224)= 0;
    mask(41, 196:200)= 0;
    mask(20:21, 246:249)= 0;
    mask(237:238, 121:123)= 0;
    mask(74:77, 91:94) = 0;
    mask(151:152, 203:207) = 0;
    mask(121:122, 209:212) = 0;
    mask(92:93, 46:48) = 0;
    mask(24, 47:49) = 0;
    mask(6:7, 109:111) = 0;
    mask(203, 43:46) = 0;
    mask(179:181, 26:28) = 0;
    mask(29:30, 199) = 0;
    mask(31:33, 225:227) = 0;
    mask(60:68, 248:251) = 0;
    mask(203, 236:237) = 0;
    mask(208, 249) = 0;
elseif strcmp(data.header.timestamp,'2458351.2088543')
    mask(156:165, 70:79)= 0;
    mask(164:169, 96:101)= 0;
    mask(32, 55:65)= 0;
    mask(48:50, 187:189)= 0;
    mask(147:148, 161:163)= 0;
    mask(166:167, 146:152)= 0;
    mask(164:166, 166:171)= 0;
    mask(32, 21:24)= 0;
    mask(36:38, 215:220)= 0;
    mask(42:43, 201:205)= 0;
    mask(236:237, 132:144)= 0;
    mask(195:198, 119:121) = 0;
    mask(202:204, 129:130) = 0;
    mask(224:226, 124:125) = 0;
    mask(241:242, 103:108) = 0;
    mask(238:240 ,110:112) = 0;
    mask(219:220, 178:183) = 0;
    mask(224:226, 204:210) = 0;
    mask(186:187, 14:17) = 0;
    mask(111:115, 206:108) = 0;
    mask(110:111, 209:210) = 0;
    mask(116:117, 211:212) = 0;
    mask(99, 212:214) = 0;
    mask(72:73, 196:198) = 0;
    mask(48:51, 64:66) = 0;
    mask(30:31, 81:87) = 0;
    mask(40:42, 19:22) = 0;
    mask(118:120, 48:50) = 0;
    mask(120:121, 52:52) = 0; 
    mask(125:126, 48:49) = 0;
    mask(127:128, 23:26) = 0;
    mask(82, 7:8) = 0;
    mask(88:92, 14:16) = 0;
    mask(7:8, 232:234) = 0;
elseif strcmp(data.header.timestamp,'2458351.2092015')
    mask(103:106, 6:19) = 0;
    mask(157:164, 70:77) = 0;
    mask(166:168, 99:101) = 0;
    mask(163:164, 103:105) = 0;
    mask(179:180, 42:46) = 0;
    mask(216:221, 109:113) = 0;
    mask(218:219, 195:201) = 0;
    mask(211:213, 192:194) = 0;
    mask(172:175, 185:189) = 0;
    mask(92:95, 166:171) = 0;
    mask(89:90, 169:179) = 0;
    mask(85, 169:171) = 0;
    mask(13:14, 49:51) = 0;
    mask(11:12, 53:59) = 0;
elseif strcmp(data.header.timestamp,'2458351.2095488')
    mask(157:165, 74:77) = 0;
    mask(165:168, 99:101) = 0;
    mask(131:133, 40:50) = 0;
    mask(71, 166:181) = 0;
    mask(62:63, 178:181) = 0;
    mask(28:29, 114:115) = 0;
    mask(153:154, 208:210) = 0;
    mask(147:148, 201:202) = 0;
    mask(149:150, 214:222) = 0;
    mask(216:217, 234:235) = 0;
    mask(224, 239:241) = 0;
    mask(166:167, 242:244)=0;
    mask(161:162, 249:251)=0;
elseif strcmp(data.header.timestamp,'2458351.2098960')
    mask(157:164, 74:78) = 0;
    mask(165:168, 98:101) = 0;
    mask(147:150, 167:177) = 0;
    mask(148:149, 181:183) = 0;
    mask(187:188, 55:63) = 0;
    mask(196:198, 120:121) = 0;
    mask(198:199, 122:127) = 0;
    mask(195:196, 127:130) = 0;
    mask(201:204, 129:130) = 0;
    mask(207:212, 114:118) = 0;
    mask(243:245, 204:213) = 0;
    mask(18:21, 96:98) = 0;
    mask(14:15, 91:95) = 0;
    mask(18:19, 107:109) = 0;
    mask(40:41, 116:117) = 0;
    mask(46:48, 129: 131) = 0;
    mask(68, 228:229) = 0;
    mask(53:54, 236:237) = 0;
    mask(57:58, 215:217) = 0;
    mask(126:127, 205:209) =0 ;
    mask(33:34, 88:92)=0;
    mask(223:224, 8:11)=0;
elseif strcmp(data.header.timestamp,'2458351.2102432')
    mask(156:164, 74:78)=0;
    mask(165:169, 98:101)=0;
    mask(130:131, 85:89)=0;
    mask(132:133, 78:80)=0;
    mask(220:221, 124:126)=0;
    mask(205:206, 50:52)=0;
    mask(200:207, 61:64)=0;
    mask(211:213, 65:70)=0;
    mask(217:218, 76:79)=0;
    mask(146:147, 227:230)=0;
    mask(101:102, 190:191)=0;
    mask(112:116, 191:192)=0;
    mask(119, 221:225)=0;
    mask(111:114, 205:209)=0;
    mask(123, 202:207)=0;
    mask(86, 246:248)=0;
    mask(91, 176:180)=0;
    mask(68:69, 214:217)=0;
    mask(68, 228:229)=0;
    mask(73, 155:161)=0;
    mask(67, 152:154)=0;
    mask(60:61, 139:141)=0;
    mask(8, 112:117)=0;
    mask(32, 23:24)=0;
    mask(17:18, 31:32)=0;
    mask(100, 34:36)=0;
    mask(135, 54:57)=0;
    mask(145:147, 26:29)=0;
    mask(163:164, 31:33)=0;
    mask(184:189, 29:31)=0;
elseif strcmp(data.header.timestamp,'2458351.2105904')
    mask(157:164, 74:82)=0;
    mask(165:168, 98:101)=0;
    mask(189, 124:126)=0;
    mask(184:186, 122:123)=0;
    mask(195:198, 120:121)=0;
    mask(201:204, 129:130)=0;
    mask(224:229, 120:124)=0;
    mask(217:218, 70:73)=0;
    mask(233:236, 181:183)=0;
    mask(177:179, 163:172)=0;
    mask(216:217, 234:237)=0;
    mask(198:204, 249:251)=0;
    mask(217:218, 70:73)=0;
    mask(160, 20:21)=0;
    mask(146:148, 27:28)=0;
    mask(134:135, 29:30)=0;
    mask(95:96, 48:50)=0;
    mask(90:92, 154:156)=0;
    mask(126:128, 219:226)=0;
    mask(109:110, 187:190)=0;
    mask(106:107, 187:191)=0;
    mask(106:107, 206:209)=0;
    mask(60:62, 196:202)=0;
    mask(71:72, 205:210)=0;
    mask(45:46, 215:217)=0;
    mask(14:15, 107:119)=0;
    mask(19:20, 136:138)=0;
elseif strcmp(data.header.timestamp,'2458351.2109377')
    mask(157:164, 74:78)=0;
    mask(164:168, 99:101)=0;
    mask(200:202, 19:21)=0;
    mask(223:224, 19:20)=0;
    mask(221:226, 52:55)=0;
    mask(200:205, 219:224)=0;
    mask(169:171, 188:192)=0;
    mask(87:92, 240:246)=0;
    mask(53:55, 39:48)=0;
    mask(75:76, 92:93)=0;
    mask(21:22, 58:60)=0;
    mask(250:251, 213:222)=0;
elseif strcmp(data.header.timestamp,'2458351.2112849')
    mask(157:164, 74:78)=0;
    mask(164:168, 99:101)=0;
    mask(17:21, 245:251)=0;
    mask(112:113, 205:207)=0;
    mask(196:198, 120:121)=0;
    mask(75:76, 92:93)=0;
% elseif strcmp(data.header.timestamp,'2458351.2116321')
%     mask();
else
    %otherwise check for interactive file
    fileout = [paths.mandir,'manmask_interactive/',data.header.timestamp,'_inter.mat']; %file name to use for interactive
    fileout_edit = [paths.mandir,'manmask_interactive/',data.header.timestamp,'_inter-e.mat']; %set as a file name
    if isfile(fileout)
        fileoutHolder = load(fileout); %load in holder (helps parfor be safe in future)
        mask = fileoutHolder.mask'; %only need this, but it also holds the direct coordinates in mask_xy
    else
        % code to help manually mask, either 'man' by typing in coordinates or
        % 'inter' which automagically does the legwork for you
        if( strcmp(manmask_mode,'man') )
            disp('\nIn debugger mode. In manual mask mode. Add a manual mask entry for the current file. Rerun code to use new manual mask.')
            ST = dbstack; dbstop('in', ST(1).file, 'at', num2str(ST(1).line+1)); %drop into debugger here
        elseif( strcmp(manmask_mode,'inter') )
            if isfile(fileout_edit) %file name if edit req'd
                %enables editing of an existing file by editing the existing file to have a -e
                fileoutHolder = load(fileout_edit); %load in holder (helps parfor be safe in future)
                mask = fileoutHolder.mask; %get mask to edit
                mask_xy = fileoutHolder.mask_xy; %need mask_xy as well as mask if editing
                mask_xy = unique(mask_xy,'rows'); %prev versions of REKT MODE accidentally allowed for duplicates to possibly occur, prevents them
            end
            %Plot side-by-side masked data without and with clip mask with
            %optional mouse movement returns value and interactive masking
            ax1_data = data.data.*(~data.mask.onemask); %with clip mask
            ax2_data = data.data.*(~data.mask.starmask & ~data.mask.linemask & ~data.mask.statmask & ~data.mask.ghostmask & ~data.mask.manmask); %without clip mask
    
            %Prep figure
            if isempty(groot().Children)
                %no figures at all, make it
                h = figure();
                clf;
                if( strcmp(manmask_mode_interOrientation,'vert') )
                    set(h,'Position',[1150 200 768 1024]); %set fig size [leftOffset bottomOffset horizontalFigSize verticalFigSize]
                elseif( strcmp(manmask_mode_interOrientation,'horiz') )
                    set(h,'Position',[200 300 1024 468]); %set fig size [leftOffset bottomOffset horizontalFigSize verticalFigSize]
                end
                h.Name = 'MaskMode'; %custom name to discern if the figure is made
                            %dynamic clicking stuff
                fprintf(['\nEntering interactive manual mask mode.\n',...
                    'Clicking with the LEFT mouse button will mask (or unmask if already masked) the selected pixel.\n',...
                    'Clicking with the RIGHT mouse button will allow you to zoom. Click the RIGHT mouse button again when done zooming to re-enter interactive mode.\n',...
                    'Clicking with the MIDDLE mouse button OR the H button will reset the zoom.\n',...
                    'Hit Q key while in the figure OR close the window to exit interactive mode and save the changes. Edit the changes by going to the file and adding -e to the file name.\n',...
                    'R activates Rectangle Mode, F resets the 1st corner in Rectangle Mode.\n',...
                    'For REKT MODE Z undoes the last Rectangle Mode selection. Y redoes it (in case undo was in error).\n',...
                    'G highlights outliers (>10, <-10) for Clip Mask Image on both images. B highlights outliers (>10, <-10) for Non-Clipped Image on both images. Usable together.',...
                    'K twice resets the ENTIRE mask (if big mistake occurs). V shows the min and max for each respective plot. Pressing again turns it off.',...
                    'This will only be printed once but can be viewed within nh_make_manmask.m script.']);
            else
                figz = groot().Children; %get figures open
                figz_right = false; %holder for the right figure "#"
                for jj = 1:length(figz)
                    if strcmp(figz(jj).Name,'MaskMode')
                        figz_right = jj; %record the right figure "#"
                    end
                end
                if figz_right == false
                    %no correct figure found, make it
                    h = figure();
                    clf; %clear it out for new stuff
                    if( strcmp(manmask_mode_interOrientation,'vert') )
                        set(h,'Position',[1150 200 768 1024]); %set fig size [leftOffset bottomOffset horizontalFigSize verticalFigSize]
                    elseif( strcmp(manmask_mode_interOrientation,'horiz') )
                        set(h,'Position',[200 300 1024 468]); %set fig size [leftOffset bottomOffset horizontalFigSize verticalFigSize]
                    end                    
                    h.Name = 'MaskMode'; %custom name to discern if the figure is made
                    fprintf(['\nEntering interactive manual mask mode.\n',...
                        'Clicking with the LEFT mouse button will mask (or unmask if already masked) the selected pixel.\n',...
                        'Clicking with the RIGHT mouse button will allow you to zoom. Click the RIGHT mouse button again when done zooming to re-enter interactive mode.\n',...
                        'Clicking with the MIDDLE mouse button OR the H button will reset the zoom.\n',...
                        'Hit Q key while in the figure OR close the window to exit interactive mode and save the changes. Edit the changes by going to the file and adding -e to the file name.\n',...
                        'R activates Rectangle Mode, F resets the 1st corner in Rectangle Mode.\n',...
                        'For REKT MODE Z undoes the last Rectangle Mode selection. Y redoes it (in case undo was in error).\n',...
                        'G highlights outliers (>10, <-10) for Clip Mask Image on both images. B highlights outliers (>10, <-10) for Non-Clipped Image on both images. Usable together.',...
                        'K twice resets the ENTIRE mask (if big mistake occurs). V shows the min and max for each respective plot. Pressing again turns it off.',...
                        'This will only be printed once but can be viewed within nh_make_manmask.m script.']);
                else
                    %figure found, use it
                    h = figz(jj);
                    clf; %clear it out for new stuff
                    fprintf('\nEntering interactive manual mask mode.\n');
                end
            end
            if( strcmp(manmask_mode_interOrientation,'vert') )
                tiledlayout(h,2,1, 'Padding', 'compact', 'TileSpacing', 'compact');
            elseif( strcmp(manmask_mode_interOrientation,'horiz') )
                tiledlayout(h,1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
            end
    
    %         ax1 = subplot(2,1,1);
            ax1 = nexttile;
            ax1_im = imagesc(ax1_data);
            %set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            % set(h,'visible','off');
            a = colorbar;
            a.Label.String = 'Intensity [DN]';
            pbaspect([1 1 1]);
            xlabel('LORRI X Pixels');
            ylabel('LORRI Y Pixels');
            caxis([-10,10]);
            set(gca,'YDir','normal');
            th1 = title(''); %place holder
            titlePos = get(th1, 'position');
            titlePosNew = titlePos + [-30, 0, 0];
            set(th1,'Position',titlePosNew);
            set(th1, 'units', 'normal');
            
    %         ax2 = subplot(2,1,2);
            ax2 = nexttile;
            ax2_im = imagesc(ax2_data);
            %set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            a = colorbar;
            a.Label.String = 'Intensity [DN]';
            pbaspect([1 1 1]);
            xlabel('LORRI X Pixels');
            ylabel('LORRI Y Pixels');
            caxis([-10,10]);
            %         title(sprintf('%s',data.header.rawfile));
            % grid minor;
            % title(sprintf('Clip-masking > %.2f + %.0f*%.2f = %.2f',clipmean,nsig,clipstd,(clipmean+nsig*clipstd)));
            % title(sprintf('Field: %d',data.header.timestamp));
            set(gca,'YDir','normal');
            th2 = title(''); %place holder
            titlePos = get(th2, 'position');
            titlePosNew = titlePos + [-30, 0, 0];
            set(th2,'Position',titlePosNew);
            set(th2, 'units', 'normal');
            
            %link axes for syncronized zoom
            linkaxes([ax1, ax2],'xy');
%             axImArray = [ax1_im, ax2_im]; %used for background mask updates (SPEED!)
            hold on;
    
            fprintf(['\nWorking on file: ',data.header.timestamp,'\n']);
            
            %prep needed vars
            if exist('mask_xy','var') == false
                mask_xy = [];
            end
            set(ax1,'color',[1 1 1]); %sets background to white so masked pixels look white
            set(ax2,'color',[1 1 1]);
            set(ax1_im, 'AlphaData', mask'); %prime the mask
            set(ax2_im, 'AlphaData', mask');
            mask_size = size(mask); %get size
            FLG_rectangleMode = false; %prime rectangle mode
            FLG_rectangleMode_c1 = false; %prime chevron 1 for rectangle mode
            rectangleMode_lastMask_undoneFLG = false; %prime rectangle mode undo
            FLG_resetL1 = false; %prime reset mode
            FLG_vMode = false; %prime reporting mode
            FLG_cMode = false; %prime caxis limit adjuster
            [mask_meshgrid_X,mask_meshgrid_Y] = meshgrid(1:1:256,1:1:256); %get a meshgrid to get X/Y indexes for what find returns
            FLG_gMode = false; %prime caxis exceeder highlighter
            FLG_bMode = false; %prime caxis exceeder highlighter

            %prime by entering zoom 1st (b/c its first thing you wanna do)
            try
                rightClick = true;
                title(ax1,'ZOOM MODE ACTIVE'); %start off in zoom mode bc makes sense
                title(ax2,'ZOOM MODE ACTIVE');
                zoom on; %make sure zoom is on
                while rightClick %loop to catch only right clicks
                    clickOrButton = waitforbuttonpress; %wait for button/click to occur
                    if clickOrButton == 0
                        if( strcmp(h.SelectionType,'alt') )
                            rightClick = false; %exit out if right clicked
                            zoom off; %make sure zoom is off
                        elseif( strcmp(h.SelectionType,'extend') )
                            zoom out;
                        end
                    else
                        %otherwise it was a button press
                        button = h.CurrentCharacter; %clicks don't register as a button, seem to be '' usually
                        if isequal(button, 'h')
                            zoom out;
                        elseif isequal(button, 'g')
                            %G mode - highlights outliers for axis 1
                            if FLG_gMode == false
                                ax1_circle = find((ax1_data > 10) | (ax1_data < -10)); %find returns flattened indexes
                                mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                hold([ax1 ax2],'on');
                                ax1_circleG_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                ax2_circleG_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                FLG_gMode = true;
                            else
                                delete(ax1_circleG_im); %remove it
                                delete(ax2_circleG_im);
                                FLG_gMode = false;
                            end
                        elseif isequal(button, 'b')
                            %B mode - highlights outliers for axis 2
                            if FLG_bMode == false
                                ax1_circle = find((ax2_data > 10) | (ax2_data < -10)); %find returns flattened indexes
                                mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                hold([ax1 ax2],'on');
                                ax1_circleB_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                ax2_circleB_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                FLG_bMode = true;
                            else
                                delete(ax1_circleB_im); %remove it
                                delete(ax2_circleB_im);
                                FLG_bMode = false;
                            end
                        elseif isequal(button, 'q')
                            % exit while loop if escape is pressed
                            disp(['Q pressed, exiting interactive masking and saving data as:\n',fileout]);
                            % save mask & mask_x/mask_y
                            if not(isfolder([paths.mandir,'manmask_interactive/']))
                                mkdir([paths.mandir,'manmask_interactive/'])
                            end
                            save(fileout,'mask','mask_xy');
        
                            break %actual exit
                        end
                    end
                end
                title(ax1,'Masking >.>'); %after leaving zoom mode back to default masking
                title(ax2,'Masking <.<');
            catch ME
                %should only activate if window is closed
                if strcmp(ME.identifier,'MATLAB:UI:CancelWaitForKeyOrButtonPress')
                    disp(['Window closed, exiting interactive masking and saving data as:\n',fileout]);
                    % save mask & mask_x/mask_y
                    if not(isfolder([paths.mandir,'manmask_interactive/']))
                        mkdir([paths.mandir,'manmask_interactive/'])
                    end
                    save(fileout,'mask','mask_xy');
                else
                    disp('RANDOM ERROR MESSAGE OCCURED:');
                    disp(ME.identifier);
                    disp(ME.message);
                    disp('dropping into debug mode');
                    ST = dbstack; dbstop('in', ST(1).file, 'at', num2str(ST(1).line+1)); %drop into debugger here
                    disp(''); %catch for debug mode b/c it catches on next line
                end
            end
            %time for interactive while loop
            while ishandle(h)
                try
                    clickOrButton = waitforbuttonpress; %wait for button/click to occur
                    if clickOrButton == 0
                        if strcmp(h.SelectionType,'normal')
                            %use point only if clicked
                            zoom off; %ensure zoom is off always
                            FLG_goodClick = false; %prime the flag
                            xy_clickAx1 = get(ax1, 'CurrentPoint');
                            xy_clickAx2 = get(ax2, 'CurrentPoint');
                            ax_xLim = ax1.XLim; %ax1 and ax2 are linked so only need to check ax1
                            ax_yLim = ax1.YLim; %ax1 and ax2 are linked so only need to check ax1
                            % if( (xy_clickAx1(1,3) == 1) && (xy_clickAx2(1,3) > 1) ) %useful thing not applicable for horizontal plots
                            if( (xy_clickAx1(1,1) >= ax_xLim(1)) && (xy_clickAx1(1,1) <= ax_xLim(2)) &&...
                                (xy_clickAx1(1,2) >= ax_yLim(1)) && (xy_clickAx1(1,2) <= ax_yLim(2)) )
                                %click inside ax1
                                x_click = round(xy_clickAx1(1,1)); %round them, the pixel 89 is shown visually from 88.5 to 89.49_ (1 is 0.5 to 1.49_)
                                y_click = round(xy_clickAx1(1,2));
                                if( (x_click <= mask_size(1)) && (x_click >= 1)  && (y_click <= mask_size(2)) && (y_click >= 1) )
                                    FLG_goodClick = true; %it's a good click if within plot bounds
                                end
                            % elseif( (xy_clickAx1(1,3) > 1) && (xy_clickAx2(1,3) == 1) )
                            elseif( (xy_clickAx2(1,1) >= ax_xLim(1)) && (xy_clickAx2(1,1) <= ax_xLim(2)) &&...
                                    (xy_clickAx2(1,2) >= ax_yLim(1)) && (xy_clickAx2(1,2) <= ax_yLim(2)) )
                                %click inside ax2
                                x_click = round(xy_clickAx2(1,1)); %round them, the pixel 89 is shown visually from 88.5 to 89.49_ (1 is 0.5 to 1.49_)
                                y_click = round(xy_clickAx2(1,2));
                                if( (x_click <= mask_size(1)) && (x_click >= 1)  && (y_click <= mask_size(2)) && (y_click >= 1) )
                                    FLG_goodClick = true; %it's a good click if within plot bounds
                                end
                            end
                            if( FLG_goodClick == true )
                                %regular click mode
                                if mask(x_click,y_click)
                                    mask(x_click,y_click) = 0; %set to 0
                                    mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                    %disp(['+Masked (',num2str(x_click),' , ',num2str(y_click),')']);
                                else
                                    mask(x_click,y_click) = 1; %reset to 1
                                    mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                    %disp(['-Unmasked (',num2str(x_click),' , ',num2str(y_click),')']);
                                end
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                                %check for rectangular selection mode (it
                                %uses the code above to show selected
                                %pixels!)
                                if FLG_rectangleMode == true
                                    if( all(FLG_rectangleMode_c1 == false) )
                                        title(ax1,'REKT MODE: CHEVRON 1 LOCKED');
                                        FLG_rectangleMode_c1 = [x_click; y_click]; %record coordinates for corner 1
                                    else
                                        % FLG_rectangleMode_c1 = [x_click; y_click]; %record coordinates for corner 2 (not needed)
                                        rectangleMode_extentX_min = min([FLG_rectangleMode_c1(1),x_click]); %get the extents
                                        rectangleMode_extentX_max = max([FLG_rectangleMode_c1(1),x_click]); %get the extents
                                        rectangleMode_extentY_min = min([FLG_rectangleMode_c1(2),y_click]); %get the extents
                                        rectangleMode_extentY_max = max([FLG_rectangleMode_c1(2),y_click]); %get the extents
                                        
                                        if( (mask(FLG_rectangleMode_c1(1),FLG_rectangleMode_c1(2)) == 1) && (mask(x_click,y_click) == 1) )
                                            %if both corners are set to 1, 1
                                            %the array
                                            rectangleMode_lastMask = zeros((rectangleMode_extentX_max-rectangleMode_extentX_min)*(rectangleMode_extentY_max-rectangleMode_extentY_min),2); %preallocate fully
                                            cntr = 1; %prep cntr
                                            mask(rectangleMode_extentX_min:rectangleMode_extentX_max,rectangleMode_extentY_min:rectangleMode_extentY_max) = 1; %set all to 1 (special!)
                                            for i = rectangleMode_extentX_min:rectangleMode_extentX_max
                                                for j = rectangleMode_extentY_min:rectangleMode_extentY_max
                                                    mask_xy(all(ismember(mask_xy, [i; j]),2),:) = []; %delete coordinates
                                                    rectangleMode_lastMask(cntr,:) = [i; j]; %record coordinates
                                                    cntr = cntr + 1; %increment
                                                end
                                            end
                                        else
                                            %otherwise 0 the array - mostly
                                            %what is intended
                                            rectangleMode_lastMask = zeros(2,2); %can't preallocate FULLY b/c can overwrite already white-spaced area that we need to ignore
                                            rectangleMode_lastMask(1,:) = FLG_rectangleMode_c1; %ensure these are there (will remove later if not supposed to be)
                                            rectangleMode_lastMask(2,:) = [x_click; y_click];
                                            for i = rectangleMode_extentX_min:rectangleMode_extentX_max
                                                for j = rectangleMode_extentY_min:rectangleMode_extentY_max
                                                    if( ~any(all(ismember(mask_xy, [i; j]),2)) )
                                                        mask_xy(end+1,:) = [i; j]; %record coordinates
                                                        rectangleMode_lastMask(end+1,:) = [i; j]; %record coordinates
                                                    end
                                                end
                                            end
                                            if( mask(FLG_rectangleMode_c1(1),FLG_rectangleMode_c1(2)) ) %remove if the pixel flipped (e.g. should stay 0)
                                                rectangleMode_lastMask(all(ismember(rectangleMode_lastMask, FLG_rectangleMode_c1),2),:) = []; %delete coordinates
                                            end
                                            if( mask(x_click,y_click) ) %remove if the pixel flipped (e.g. should stay 0)
                                                rectangleMode_lastMask(all(ismember(rectangleMode_lastMask, [x_click; y_click]),2),:) = []; %delete coordinates
                                            end
                                            mask(rectangleMode_extentX_min:rectangleMode_extentX_max,rectangleMode_extentY_min:rectangleMode_extentY_max) = 0; %set all to 0
                                        end
                                        %update plot dynamically to see clicked pixels
                                        set(ax1_im, 'AlphaData', mask');
                                        set(ax2_im, 'AlphaData', mask');
                                        FLG_rectangleMode_c1 = false; %reset the flag
                                        title(ax1,'REKT MODE ACTIVE'); %reset title
                                        rectangleMode_lastMask_undoneFLG = false; %prime this b/c overwritten
                                    end
                                end
                            end
                        elseif strcmp(h.SelectionType,'alt')
                            if all(FLG_rectangleMode_c1 ~= false)
                                %undo last click (which was setting RECT MODE corner 1)
                                if mask(x_click,y_click)
                                    mask(x_click,y_click) = 0; %set to 0
                                    mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                else
                                    mask(x_click,y_click) = 1; %reset to 1
                                    mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                end
                                title(ax1,'REKT MODE ACTIVE');
                                FLG_rectangleMode_c1 = false; %set to false
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                            end
                            rightClick = true;
                            rightClick_quit = false; %for breaking out of 2 while loops
                            ax1_prevTitle = ax1.Title.String;
                            ax2_prevTitle = ax2.Title.String;
                            title(ax1,'ZOOM MODE ACTIVE');
                            title(ax2,'ZOOM MODE ACTIVE');
                            zoom on; %make sure zoom is on
                            while rightClick %loop to catch only right clicks
                                clickOrButton = waitforbuttonpress; %wait for button/click to occur
                                if clickOrButton == 0
                                    if( strcmp(h.SelectionType,'alt') )
                                        rightClick = false; %exit out if right clicked
                                        zoom off; %make sure zoom is off
                                    elseif( strcmp(h.SelectionType,'extend') )
                                        zoom out;
                                    end
                                else
                                    %otherwise it was a button press
                                    button = h.CurrentCharacter; %clicks don't register as a button, seem to be '' usually
                                    if isequal(button, 'h')
                                        zoom out;
                                    elseif isequal(button, 'g')
                                        %G mode - highlights outliers for axis 1
                                        if FLG_gMode == false
                                            ax1_circle = find((ax1_data > 10) | (ax1_data < -10)); %find returns flattened indexes
                                            mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                            mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                            hold([ax1 ax2],'on');
                                            ax1_circleG_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                            ax2_circleG_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                            FLG_gMode = true;
                                        else
                                            delete(ax1_circleG_im); %remove it
                                            delete(ax2_circleG_im);
                                            FLG_gMode = false;
                                        end
                                    elseif isequal(button, 'b')
                                        %B mode - highlights outliers for axis 2
                                        if FLG_bMode == false
                                            ax1_circle = find((ax2_data > 10) | (ax2_data < -10)); %find returns flattened indexes
                                            mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                            mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                            hold([ax1 ax2],'on');
                                            ax1_circleB_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                            ax2_circleB_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                            FLG_bMode = true;
                                        else
                                            delete(ax1_circleB_im); %remove it
                                            delete(ax2_circleB_im);
                                            FLG_bMode = false;
                                        end
                                    elseif isequal(button, 'q')
                                        % exit while loop if escape is pressed
                                        disp(['Q pressed, exiting interactive masking and saving data as:\n',fileout]);
                                        % save mask & mask_x/mask_y
                                        if not(isfolder([paths.mandir,'manmask_interactive/']))
                                            mkdir([paths.mandir,'manmask_interactive/'])
                                        end
                                        save(fileout,'mask','mask_xy');
                                        
                                        rightClick_quit = true; %flip it
                                        break %actual exit (but not really)
                                    end
                                end
                            end
                            if( rightClick_quit == true )
                                break %actual exit
                            end
                            title(ax1,ax1_prevTitle);
                            title(ax2,ax2_prevTitle);
                        elseif strcmp(h.SelectionType,'extend')
                            zoom out;
                            if all(FLG_rectangleMode_c1 ~= false)
                                %undo last click (which was setting RECT MODE corner 1)
                                if mask(x_click,y_click)
                                    mask(x_click,y_click) = 0; %set to 0
                                    mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                else
                                    mask(x_click,y_click) = 1; %reset to 1
                                    mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                end
                                title(ax1,'REKT MODE ACTIVE');
                                FLG_rectangleMode_c1 = false; %set to false
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                            end
                            rightClick = true;
                            rightClick_quit = false;
                            ax1_prevTitle = ax1.Title.String;
                            ax2_prevTitle = ax2.Title.String;
                            title(ax1,'ZOOM MODE ACTIVE');
                            title(ax2,'ZOOM MODE ACTIVE');
                            zoom on; %make sure zoom is on
                            while rightClick %loop to catch only right clicks
                                clickOrButton = waitforbuttonpress; %wait for button/click to occur
                                if clickOrButton == 0
                                    if( strcmp(h.SelectionType,'alt') )
                                        rightClick = false; %exit out if right clicked
                                        zoom off; %make sure zoom is off
                                    elseif( strcmp(h.SelectionType,'extend') )
                                        zoom out;
                                    end
                                else
                                    %otherwise it was a button press
                                    button = h.CurrentCharacter; %clicks don't register as a button, seem to be '' usually
                                    if isequal(button, 'h')
                                        zoom out;
                                    elseif isequal(button, 'g')
                                        %G mode - highlights outliers for axis 1
                                        if FLG_gMode == false
                                            ax1_circle = find((ax1_data > 10) | (ax1_data < -10)); %find returns flattened indexes
                                            mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                            mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                            hold([ax1 ax2],'on');
                                            ax1_circleG_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                            ax2_circleG_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                            FLG_gMode = true;
                                        else
                                            delete(ax1_circleG_im); %remove it
                                            delete(ax2_circleG_im);
                                            FLG_gMode = false;
                                        end
                                    elseif isequal(button, 'b')
                                        %B mode - highlights outliers for axis 2
                                        if FLG_bMode == false
                                            ax1_circle = find((ax2_data > 10) | (ax2_data < -10)); %find returns flattened indexes
                                            mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                            mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                            hold([ax1 ax2],'on');
                                            ax1_circleB_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                            ax2_circleB_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                            FLG_bMode = true;
                                        else
                                            delete(ax1_circleB_im); %remove it
                                            delete(ax2_circleB_im);
                                            FLG_bMode = false;
                                        end
                                    elseif isequal(button, 'q')
                                        % exit while loop if escape is pressed
                                        disp(['Q pressed, exiting interactive masking and saving data as:\n',fileout]);
                                        % save mask & mask_x/mask_y
                                        if not(isfolder([paths.mandir,'manmask_interactive/']))
                                            mkdir([paths.mandir,'manmask_interactive/'])
                                        end
                                        save(fileout,'mask','mask_xy');
                                        
                                        rightClick_quit = true; %flip it
                                        break %actual exit (but not really)
                                    end
                                end
                            end
                            if( rightClick_quit == true )
                                break %actual exit
                            end
                            title(ax1,ax1_prevTitle);
                            title(ax2,ax2_prevTitle);
                        end
                        FLG_resetL1 = false; %keep off
                    else
                        %otherwise it was a button press
                        button = h.CurrentCharacter; %clicks don't register as a button, seem to be '' usually
                        if isequal(button, 'r')
                            %activates or deactivates rectangle mode
                            FLG_rectangleMode = ~FLG_rectangleMode; %flip the flag
                            if( FLG_rectangleMode == true )
                                title(ax1,'REKT MODE ACTIVE');
                                title(ax2,'REKT MODE ACTIVE');
                                FLG_rectangleMode_c1 = false; % prime corner 1 and corner 2 flags
                            else
                                title(ax1,'Masking >.>');
                                title(ax2,'Masking <.<');
                                if all(FLG_rectangleMode_c1 ~= false)
                                    %undo last click (which was setting RECT MODE corner 1)
                                    if mask(x_click,y_click)
                                        mask(x_click,y_click) = 0; %set to 0
                                        mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                    else
                                        mask(x_click,y_click) = 1; %reset to 1
                                        mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                    end
                                    title(ax1,'REKT MODE ACTIVE');
                                    FLG_rectangleMode_c1 = false; %set to false
                                    %update plot dynamically to see clicked pixels
                                    set(ax1_im, 'AlphaData', mask');
                                    set(ax2_im, 'AlphaData', mask');
                                end
                            end
                            FLG_resetL1 = false; %keep off
                        elseif isequal(button, 'f')
                            if all(FLG_rectangleMode_c1 ~= false)
                                disp('F-REKT MODE RESET')
                                title(ax1,'REKT MODE ACTIVE');
                                FLG_rectangleMode_c1 = false; % prime these
                                %undo last click (which was setting RECT MODE corner 1)
                                if mask(x_click,y_click)
                                    mask(x_click,y_click) = 0; %set to 0
                                    mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                else
                                    mask(x_click,y_click) = 1; %reset to 1
                                    mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                end
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                            end
                        elseif isequal(button, 'z')
                            if (FLG_rectangleMode == true) && exist('rectangleMode_lastMask','var') && (rectangleMode_lastMask_undoneFLG == false) %only run if rectangleMode_lastMask exists
                                rectangleMode_lastMask_entries = size(rectangleMode_lastMask);
                                rectangleMode_lastMask_entries = rectangleMode_lastMask_entries(1); %get the entries
                                for i = 1:rectangleMode_lastMask_entries
                                    x_click = rectangleMode_lastMask(i,1); %alias to copy code ez
                                    y_click = rectangleMode_lastMask(i,2);
                                    if mask(x_click,y_click)
                                        mask(x_click,y_click) = 0; %set to 0
                                        mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                    else
                                        mask(x_click,y_click) = 1; %reset to 1
                                        mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                    end
                                end
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                                rectangleMode_lastMask_undoneFLG = true; %set flag to allow redo
                                disp(['Z-REKT MODE: Undid last ',num2str(rectangleMode_lastMask_entries),' values selected with REKT MODE. You can undo this by pressing Y.'])
                            elseif (FLG_rectangleMode == true) && exist('rectangleMode_lastMask','var') && (rectangleMode_lastMask_undoneFLG == true)
                                disp('Z-REKT MODE: Already undid last values. You can press Y to redo those values.')
                            elseif (FLG_rectangleMode == true) && ~exist('rectangleMode_lastMask','var')
                                disp('Z-REKT MODE: No REKT MODE masking has occured yet. Nothing to undo!')
                            elseif (FLG_rectangleMode == false)
                                disp('Z-Must be in REKT MODE to undo last REKT MODE selection.')
                            end
                        elseif isequal(button, 'y')
                            if (FLG_rectangleMode == true) && exist('rectangleMode_lastMask','var') && (rectangleMode_lastMask_undoneFLG == true) %only run if rectangleMode_lastMask exists
                                rectangleMode_lastMask_entries = size(rectangleMode_lastMask);
                                rectangleMode_lastMask_entries = rectangleMode_lastMask_entries(1); %get the entries
                                for i = 1:rectangleMode_lastMask_entries
                                    x_click = rectangleMode_lastMask(i,1); %alias to copy code ez
                                    y_click = rectangleMode_lastMask(i,2);
                                    if mask(x_click,y_click)
                                        mask(x_click,y_click) = 0; %set to 0
                                        mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                    else
                                        mask(x_click,y_click) = 1; %reset to 1
                                        mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                    end
                                end
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                                rectangleMode_lastMask_undoneFLG = false; %set flag to allow undo
                                disp(['Y-REKT MODE: Redid last ',num2str(rectangleMode_lastMask_entries),' values selected with REKT MODE. You can re-undo this by pressing Z.'])
                            elseif (FLG_rectangleMode == true) && exist('rectangleMode_lastMask','var') && (rectangleMode_lastMask_undoneFLG == false) 
                                disp('Y-REKT MODE: Cannot redo selection because nothing has been undone (or a new REKT MODE selection overwrote the previous undo).')
                            elseif (FLG_rectangleMode == true) && ~exist('rectangleMode_lastMask','var')
                                disp('Y-REKT MODE: No REKT MODE masking has occured yet. Nothing to possibly redo after being undone!')
                            elseif (FLG_rectangleMode == false)
                                disp('Y-Must be in REKT MODE to redo last REKT MODE selection (if there is one).')
                            end
                        elseif isequal(button, 'k')
                            if( FLG_resetL1 == false )
                                disp('K-WARNING: RESET LEVEL 1 ACTIVATED, PRESS AGAIN TO RESET ENTIRE MASK');
                                FLG_resetL1 = true;
                            else
                                disp('K-WARNING: WHOLE MASK RESET ACTIVATED');
                                FLG_resetL1 = false;
                                mask_xy = []; %delete all
                                mask(:,:) = 1; %set to 1
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                            end
                        elseif isequal(button, 'h')
                            zoom out;
                            if all(FLG_rectangleMode_c1 ~= false)
                                %undo last click (which was setting RECT MODE corner 1)
                                if mask(x_click,y_click)
                                    mask(x_click,y_click) = 0; %set to 0
                                    mask_xy(end+1,:) = [x_click; y_click]; %record coordinates
                                else
                                    mask(x_click,y_click) = 1; %reset to 1
                                    mask_xy(all(ismember(mask_xy, [x_click; y_click]),2),:) = []; %delete coordinates
                                end
                                title(ax1,'REKT MODE ACTIVE');
                                FLG_rectangleMode_c1 = false; %set to false
                                %update plot dynamically to see clicked pixels
                                set(ax1_im, 'AlphaData', mask');
                                set(ax2_im, 'AlphaData', mask');
                            end
                            rightClick = true;
                            rightClick_quit = false;
                            ax1_prevTitle = ax1.Title.String;
                            ax2_prevTitle = ax2.Title.String;
                            title(ax1,'ZOOM MODE ACTIVE');
                            title(ax2,'ZOOM MODE ACTIVE');
                            zoom on; %make sure zoom is on
                            while rightClick %loop to catch only right clicks
                                clickOrButton = waitforbuttonpress; %wait for button/click to occur
                                if clickOrButton == 0
                                    if( strcmp(h.SelectionType,'alt') )
                                        rightClick = false; %exit out if right clicked
                                        zoom off; %make sure zoom is off
                                    elseif( strcmp(h.SelectionType,'extend') )
                                        zoom out;
                                    end
                                else
                                    %otherwise it was a button press
                                    button = h.CurrentCharacter; %clicks don't register as a button, seem to be '' usually
                                    if isequal(button, 'h')
                                        zoom out;
                                    elseif isequal(button, 'g')
                                        %G mode - highlights outliers for axis 1
                                        if FLG_gMode == false
                                            ax1_circle = find((ax1_data > 10) | (ax1_data < -10)); %find returns flattened indexes
                                            mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                            mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                            hold([ax1 ax2],'on');
                                            ax1_circleG_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                            ax2_circleG_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                            FLG_gMode = true;
                                        else
                                            delete(ax1_circleG_im); %remove it
                                            delete(ax2_circleG_im);
                                            FLG_gMode = false;
                                        end
                                    elseif isequal(button, 'b')
                                        %B mode - highlights outliers for axis 2
                                        if FLG_bMode == false
                                            ax1_circle = find((ax2_data > 10) | (ax2_data < -10)); %find returns flattened indexes
                                            mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                            mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                            hold([ax1 ax2],'on');
                                            ax1_circleB_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                            ax2_circleB_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                            FLG_bMode = true;
                                        else
                                            delete(ax1_circleB_im); %remove it
                                            delete(ax2_circleB_im);
                                            FLG_bMode = false;
                                        end
                                    elseif isequal(button, 'q')
                                        % exit while loop if escape is pressed
                                        disp(['Q pressed, exiting interactive masking and saving data as:\n',fileout]);
                                        % save mask & mask_x/mask_y
                                        if not(isfolder([paths.mandir,'manmask_interactive/']))
                                            mkdir([paths.mandir,'manmask_interactive/'])
                                        end
                                        save(fileout,'mask','mask_xy');
                                        
                                        rightClick_quit = true; %flip it
                                        break %actual exit (but not really)
                                    end
                                end
                            end
                            if( rightClick_quit == true )
                                break %actual exit
                            end
                            title(ax1,ax1_prevTitle);
                            title(ax2,ax2_prevTitle);
                            FLG_resetL1 = false; %keep off
                        elseif isequal(button, 'v')
                            FLG_vMode = ~FLG_vMode; %flip the flag
                            if FLG_vMode
                                ax1_prevTitle_forV = ax1.Title.String;
                                ax2_prevTitle_forV = ax2.Title.String;
                                title(ax1,[ax1_prevTitle_forV,'|MAX: ',num2str(round(max(max(ax1_data.*mask)),5)),' | ','MIN: ',num2str(round(min(min(ax1_data.*mask)),5))]);
                                title(ax2,[ax2_prevTitle_forV,'|MAX: ',num2str(round(max(max(ax2_data.*mask)),5)),' | ','MIN: ',num2str(round(min(min(ax2_data.*mask)),5))]);
                            else
                                if ~isempty(strfind(ax1.Title.String,'MAX: '))
                                    %return to original
                                    title(ax1,ax1_prevTitle_forV);
                                    title(ax2,ax2_prevTitle_forV);
                                else
                                    %otherwise something overrode it, so
                                    %just reactivate b/c person wants to
                                    %see
                                    FLG_vMode = true;
                                    ax1_prevTitle_forV = ax1.Title.String;
                                    ax2_prevTitle_forV = ax2.Title.String;
                                    title(ax1,[ax1_prevTitle_forV,'|MAX: ',num2str(round(max(max(ax1_data.*mask)),5)),' | ','MIN: ',num2str(round(min(min(ax1_data.*mask)),5))]);
                                    title(ax2,[ax2_prevTitle_forV,'|MAX: ',num2str(round(max(max(ax2_data.*mask)),5)),' | ','MIN: ',num2str(round(min(min(ax2_data.*mask)),5))]);
                                end
                            end
                            FLG_resetL1 = false; %keep off
                        elseif isequal(button, 'c')
                            %C mode - disables caxis limits to see outliers
                            %really obviously (hopefully)
                            if FLG_cMode == false
                                caxis(ax1, [-inf, inf])
                                caxis(ax2, [-inf, inf])
                                FLG_cMode = true;
                            else
                                caxis(ax1, [-10, 10])
                                caxis(ax2, [-10, 10])
                                FLG_cMode = false;
                            end
                        elseif isequal(button, 'g')
                            %G mode - highlights outliers for axis 1
                            if FLG_gMode == false
                                ax1_circle = find((ax1_data > 10) | (ax1_data < -10)); %find returns flattened indexes
                                mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                hold([ax1 ax2],'on');
                                ax1_circleG_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                ax2_circleG_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'r','x','LineWidth',3); %X
                                FLG_gMode = true;
                            else
                                delete(ax1_circleG_im); %remove it
                                delete(ax2_circleG_im);
                                FLG_gMode = false;
                            end
                        elseif isequal(button, 'b')
                            %B mode - highlights outliers for axis 2
                            if FLG_bMode == false
                                ax1_circle = find((ax2_data > 10) | (ax2_data < -10)); %find returns flattened indexes
                                mask_circle_X = mask_meshgrid_X(ax1_circle); %get X indexes
                                mask_circle_Y = mask_meshgrid_Y(ax1_circle); %get X indexes
                                hold([ax1 ax2],'on');
                                ax1_circleB_im = scatter(ax1,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                ax2_circleB_im = scatter(ax2,mask_circle_X,mask_circle_Y,55,'m','+','LineWidth',3); %+
                                FLG_bMode = true;
                            else
                                delete(ax1_circleB_im); %remove it
                                delete(ax2_circleB_im);
                                FLG_bMode = false;
                            end
                        elseif isequal(button, 'q')
                            % exit while loop if escape is pressed
                            disp(['Q pressed, exiting interactive masking and saving data as:\n',fileout]);
                            % save mask & mask_x/mask_y
                            if not(isfolder([paths.mandir,'manmask_interactive/']))
                                mkdir([paths.mandir,'manmask_interactive/'])
                            end
                            save(fileout,'mask','mask_xy');
        
                            break %actual exit
                        else
                            FLG_resetL1 = false; %keep off
                        end
                    end
                    if ~isempty(strfind(ax1.Title.String,'MAX: '))
                        % Update V mode if it is on
                        title(ax1,[ax1_prevTitle_forV,'|MAX: ',num2str(round(max(max(ax1_data.*mask)),5)),' | ','MIN: ',num2str(round(min(min(ax1_data.*mask)),5))]);
                        title(ax2,[ax2_prevTitle_forV,'|MAX: ',num2str(round(max(max(ax2_data.*mask)),5)),' | ','MIN: ',num2str(round(min(min(ax2_data.*mask)),5))]);
                    end
                catch ME
                    %should only activate if window is closed
                    if strcmp(ME.identifier,'MATLAB:UI:CancelWaitForKeyOrButtonPress')
                        disp(['Window closed, exiting interactive masking and saving data as:\n',fileout]);
                        % save mask & mask_x/mask_y
                        if not(isfolder([paths.mandir,'manmask_interactive/']))
                            mkdir([paths.mandir,'manmask_interactive/'])
                        end
                        save(fileout,'mask','mask_xy');
                    else
                        disp('RANDOM ERROR MESSAGE OCCURED:');
                        disp(ME.identifier);
                        disp(ME.message);
                        disp('dropping into debug mode');
                        ST = dbstack; dbstop('in', ST(1).file, 'at', num2str(ST(1).line+1)); %drop into debugger here
                        disp(''); %catch for debug mode b/c it catches on next line
                    end
                end
            end
            if isfile(fileout_edit) %file name if edit occured
                delete(fileout_edit); %remove -e file b/c file saved successfully
            end
            % ext = '.png';
            % if not(isfolder([paths.clipmaskdir]))
            %     mkdir([paths.clipmaskdir])
            % end
            % imagename = sprintf('%s%s%s',paths.clipmaskdir,data.header.timestamp,ext);
            % imagename = sprintf('%s%s%s',paths.selectmaskdir,data.header.timestamp,ext);
            % print(h,imagename, '-dpng');
        end
        mask = mask'; %transpose for use
    end
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
