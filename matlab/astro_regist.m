%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C.Nguyen
% MATLAB routine to solve the field of images, by producing a string that
% acts as command line input for solve_field
% 11/18/2015: ver 1.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%
% nomra = initial guess for RA
% nomdec = initial guess for DEC
% solf = array with parameters of solve_field. Basically what this script
% will do is put in some initial guesses for these parameters
% inst = index for solf(). inst = 1 or 2 (they are identical)
%
% imfile = path to input image file to be performed solve_field on.
% befile = path to back end file. Back end file will point solve_field to
% the index library where solve_field can match stars to.
% otfile = path to .axy file that solve_field will create and then
% reference as it solves for fields. This file is temporary and can be
% cleaned up later.
% offile = path to result .fit file. The information we need is in the
% header of this file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function astro_regist(phase,path)
% load the RA and DEC data for each image

load (sprintf('%smatrix/%s_ra_dec.mat',path,phase));

[Nrows,Ncols] = size(RA_array);

i = 1;

for i=1:Ncols
    sprintf('Currently on file number %d',i)
    
    nomra = RA_array(i);
    nomdec = DEC_array(i);
    filename = file_name_array{i};
    
    % get the field solver parameters
    solf = get_solvefields_params();
    inst = 1;
    solf(1).defrad = '1';
    
    % make some paths
    imfile = sprintf('%s%s/selected_data/%s',path,phase,filename);
    %befile = sprintf('mybackend.cfg');
    befile = sprintf('/usr/local/astrometry/etc/astrometry.cfg')
    otfile = sprintf('reduc_astro_align');
    offile = sprintf('%sregist/selected_data/regist_%s_%s',path,phase,filename);
    
    solfcall = sprintf(...
        ['solve-field %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s '...
        ' %s %f %s %f %s %s %s %s %s %s %s'],...
        solf(inst).overwrite,solf(inst).noplots,...
        solf(inst).backend,befile,...
        solf(inst).outdir,sprintf('%sregist/',path),...
        solf(inst).outname,otfile,...
        solf(inst).fitsname,offile,...
        solf(inst).strscalelow,solf(inst).scalelow,...
        solf(inst).strscalehigh,solf(inst).scalehigh,...
        solf(inst).strscaleunits,solf(inst).scaleunits,...
        solf(inst).ra,nomra,...
        solf(inst).dec,nomdec,...
        solf(inst).radius,solf(inst).defrad,...
        solf(inst).tweak,solf(inst).tord,solf(inst).center,...
        imfile);
    
    % now the solver's been set up, do the system call to
    s = system(solfcall);
    
    % remove all of the extra files that solve-field makes - these are not
    % interesting and just take up space.
    
    
    
    rmcall = sprintf('rm %sregist/%s*',path, otfile);
    s = system(rmcall);
    
    %% now read in the header info from the solved field - we need to parse
    %% this to get out the relevant info for matlab
    
    % check to see if output file exists meaning the field was successfully
    % solved - if field was not solved skip file and move on to next image
    if isfile(offile)
        
        offile
        astro_head = fitsinfo(offile);
        
        %% pass it to a subroutine that creates the astrometry structure for
        %% this field
        
        astrometry = parse_fits_astrom_header(astro_head);
        
        % append astrometry to data structure
        data.astrometry = astrometry;
        
        
        %% save the new data
        outfile = sprintf('%sregist/%s_%s.mat',path, phase, filename);
        save(outfile,'data');
        
    else
        sprintf('Output file does not exist, skipping to next image %d.',i+1)
    end
    i = i+1
    
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function astrometry = parse_fits_astrom_header(astro_head)

header = astro_head.PrimaryData.Keywords;
nentries = length(header);

astrometry.header = header;
%
for ihead=1:nentries
    thisentry = lower(header{ihead,1});
    if ~strcmp(thisentry,'comment') & ...
            ~strcmp(thisentry,'history') & ~strcmp(thisentry,'end')
        % the information in astrometry.net header files sometimes starts
        % with '_blahblah', but setfield() requires field to start with a letter,
        % so add 'id_' to the front of the field.
        thisentry = sprintf('id_%s', thisentry);
        astrometry = setfield(astrometry,thisentry,header{ihead,2});
    end
end

[X,Y] = meshgrid([1:astrometry.id_naxis1],[1:astrometry.id_naxis2]);

X = transpose(X);
Y = transpose(Y);

[ra,dec] = pix2radec(astrometry,X(:),Y(:));
%
astrometry.ra  = (reshape(ra,astrometry.id_naxis1,astrometry.id_naxis2));
astrometry.dec = (reshape(dec,astrometry.id_naxis1,astrometry.id_naxis2));

return

