%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  function nh_pipeline.m
%%  Jun 2016, MZ
%%  For each file in the NH data directories, this program:
%%   1) Reads it in.
%%   2)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nh_pipeline()
clear all
close all

procstepflag = 1.5;

% paths = get_paths();
npaths = get_paths_new();

ndatafiles = dir(sprintf('%s*.mat',npaths.datadir));

mydate = zeros(size(ndatafiles));
mytemp = zeros(size(ndatafiles));
mymean = zeros(size(ndatafiles));
myref = zeros(size(ndatafiles));
myeng = zeros(size(ndatafiles));
myisl = zeros(size(ndatafiles));

for ifile=1:size(ndatafiles,1)
    
    disp(sprintf('On file %d of %d.',ifile,size(ndatafiles,1)));
    
    load(sprintf('%s%s',npaths.datadir,ndatafiles(ifile).name));
    
    if procstepflag > 1
        
%         data = nh_findghost(data);
        data = nh_makemask(data);
        
%        save(sprintf('%s%s',paths.datadir,ndatafiles(ifile).name),'data');
        
    end
    
    if procstepflag > 2
        
        data = nh_add_meta(data);
        
      %  save(sprintf('%s%s',paths.datadir,ndatafiles(ifile).name),'data');
        
    end
    
    if procstepflag > 3
        
        data = nh_calibrate(data);
        
       % save(sprintf('%s%s',paths.datadir,ndatafiles(ifile).name),'data');
        
    end
    
    if procstepflag > 4
        
        data = nh_calcisl(data);
        
        %save(sprintf('%s%s',paths.datadir,ndatafiles(ifile).name),'data');
        
    end
    
    if procstepflag > 5
        
        data = nh_calcdgl(data);
        
        %save(sprintf('%s%s',paths.datadir,ndatafiles(ifile).name),'data');
        
    end
    
    %     mydate(ifile) = data.header.date_jd;
    %     mytemp(ifile) = data.header.ccdtemp;
    %     mymean(ifile) = data.stats.calmean;
    %     myref(ifile) = data.ref.mean;
    %     myeng(ifile) = data.ref.engmean;
    %     %myisl(ifile) = data.isl.trimean;
    %
    %     if 0
    %     figure(1); clf
    %     imagesc(data.image.calimage)
    %     caxis([0,3000])
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	'(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_raw.png',data.header.timestamp));
    %     figure(2); clf
    %     imagesc(data.mask.mask)
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	  '(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_mask.png',data.header.timestamp));
    %     figure(3); clf
    %     imagesc(data.image.calimage.*~data.mask.mask)
    %     caxis([0,1500])
    %     colorbar
    %     title(sprintf(['%s, %5.3f, %5.3f, %5.3f\nField %d, %s, %3.1f\n'...
    % 	'(l,b) = (%6.2f, %6.2f), (elon,elat)=(%6.2f, %6.2f)'],...
    % 	data.header.timestamp,...
    % 	data.mask.maskfrac,data.stats.maskmean,data.stats.maskstd,...
    % 	data.header.fieldnum,...
    % 	data.header.target_name,data.astrometry.id_exptime,...
    % 	data.coords.galactic(1),data.coords.galactic(2),...
    % 	data.coords.ecliptic(1),data.coords.ecliptic(2)));
    %     screen2png(sprintf('plots/%s_masked.png',data.header.timestamp));
    %       end
    
    
    %dbstop
    
    
end

%figure(4);
%plot(mydate - data.header.launch_jd,mytemp,'o');

dbstop




end
