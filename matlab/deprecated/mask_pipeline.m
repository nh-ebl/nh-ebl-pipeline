function mask_pipeline(phase,path)
%
sprintf('This is phase %s',phase)

% Record JD 
sprintf('Extracting JD data...')
JD_date(phase,path);
sprintf('JD Extraction completed.')

% Find stellar flux cutoff
sprintf('Calculating stellar flux cutoff...')
cutoff_loop(phase,path);
sprintf('Stellar flux cutoff calculation completed.')

% Mask the stars
sprintf('Masking stars...')
mask(phase,path);
sprintf('Star masking completed.')

% Done
sprintf('All tasks completed. Goodbye')

end

