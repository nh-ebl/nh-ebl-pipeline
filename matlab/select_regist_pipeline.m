function select_regist_pipeline(phase1)

% path to main data folder
path = '/data/symons/nh_data_new/';

phase_array = {phase1};
size(phase_array);

for i=1:1
    phase = char(phase_array(i));
    
    print_phase(phase)
    
    
%     %construct duration data
    sprintf('Collecting exposure time data...')
%     duration_cust(phase,path);
%     texp_vs_t_ext(phase,path);
%     texp_vs_t(phase,path); % Not needed
    sprintf('Duration task completed.')
%     
%     % select data
    sprintf('Selecting useful data...')
%     selection(phase,path);
    sprintf('Selection task completed.')
    
    % construct .eps
    % sprintf('Constructing .eps images...')
    % eps_img(phase,-5,50,path);
    % sprintf('EPS construction task completed.')
    
%     % construct .png
%     %expecting files in 'masked' in phase folder
%     sprintf('Constructing .png images...')
%     %png_img(phase,-5,50,path);
%     png_img_cust(phase,-5,50,path);
%     sprintf('PNG construction task completed.')
    %retrieve header info and create png images with info annotated
    add_info(phase, path)
%     
    % extract distance information
    sprintf('Extracting trajectory information...')
%     distance(phase,path)
    % construct .txt of distance
%     t_vs_d(phase,path)
    sprintf('Distance task completed.')
    
    % register RA and dec
    sprintf('Performing astrometric registration...')
%     ra_dec(phase,path)
%     ra_dec_mat_to_dat(phase,path)
%     astro_regist(phase,path)
    sprintf('Registration task completed')
    
    % mask_pipeline(sprintf('regist'),path)
    
    % finish
    sprintf('All tasks completed. Goodbye.')
end

