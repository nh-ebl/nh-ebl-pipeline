function data= nh_dgl_manmask(data,paths)

% Display interactive image for mask determination
% figure(1); clf
% imagesc(data.dgl.ohmim_iris)
% colorbar
% pbaspect([1 1 1]);
% set(gca,'YDir','normal');
% hold on
% set(gcf, 'WindowButtonMotionFcn', @mouseMove); %overlay of mouse coords in title

mask = ones(256,256);
% data.header.timestamp

% field 5 - this field needs masking
% if strcmp(data.header.timestamp,'2457582.5835098') - first image needing
% masking
if data.header.fieldnum == 5
    skip = 0; % Default skip state
    if (str2double(data.header.timestamp) > 2456000) && (str2double(data.header.timestamp) < 2457583)
        %mask a circle of radius bigrad at [circ_x, circ_y] coords
        bigrad = 65; % radius of masking circle
        circ_x = 59; % x coordinate to center masking circle on
        circ_y = 69; % y coordinate to center masking circle on
    elseif (str2double(data.header.timestamp) > 2457583) && (str2double(data.header.timestamp) < 2457583.507)
        %mask a circle of radius bigrad at [circ_x, circ_y] coords
        bigrad = 65; % radius of masking circle
        circ_x = 105; % x coordinate to center masking circle on
        circ_y = 91; % y coordinate to center masking circle on
    elseif (str2double(data.header.timestamp) > 2457583.507)
        %mask a circle of radius bigrad at [circ_x, circ_y] coords
        bigrad = 65; % radius of masking circle
        circ_x = 88; % x coordinate to center masking circle on
        circ_y = 85; % y coordinate to center masking circle on
    else
        skip = 1; % Old data field 5 doesn't need masks
    end
    %quick n easy function to mask a circle of radius bigrad at [circ_x, circ_y] coords
%     [Xr, Yr] = meshgrid(1:1:data.astrom.imagew,1:1:data.astrom.imageh);
%     mask = sqrt((Xr-circ_x).^2+(Yr-circ_y).^2)>bigrad;
    
    if skip ~= 1
        %create submask of object (create array of 0's and 1's,
        %where the 1's represnt to location of the objects in the
        %submask). Basically a sqaure of 0's with a circle of 1's
        [rr, cc] = meshgrid(1:2*bigrad+1);
        Circle = sqrt((rr-bigrad-1).^2+(cc-bigrad-1).^2)<=bigrad;

        %combined the submask (C) and mask (Z) where xpix, ypix is
        %the center of the object
        for i = 1:(2*bigrad+1)
            xcurr = round(circ_x-bigrad-1+i);
            if xcurr < 1 || xcurr > data.astrom.imagew;
                continue;
            end
            for j = 1:(2*bigrad+1)
                ycurr = round(circ_y-bigrad-1+j);
                if ycurr < 1 || ycurr > data.astrom.imageh;
                    continue;
                end
                if Circle(i,j) == 1
                    mask(ycurr,xcurr) = 0;
                end
            end
        end
    end
    
end

% Save masked image
% h = figure(2); clf
% imagesc(data.dgl.ohmim_iris.*mask)
% colorbar
% pbaspect([1 1 1]);
% set(gca,'YDir','normal');
% hold on
% % set(gcf, 'WindowButtonMotionFcn', @mouseMove); %overlay of mouse coords in title
% ext = '.png';
% imagename = sprintf('%s%s%s',paths.dglmaskdir,data.header.timestamp,ext);
% print(h,imagename, '-dpng');

manmask = ~mask;

fileout = sprintf('%s%s.mat',paths.dglmandir,data.header.timestamp);
save(fileout,'manmask');

end