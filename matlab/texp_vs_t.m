%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C.Nguyen
% MATLAB script to export MATLAB matrix into .dat file if needed
 function texp_vs_t(phase,path)

launch_date = 2453755.291667;
au = 1.496*10^8;

filename = fopen(sprintf('%stxt/%s_texp_vs_t.txt',path,phase), 'w');
fprintf(filename, '%%t t_from_launch texp\n');

load (sprintf('%smatrix/%s_texp_t.mat',path,phase));

date = date_array;
texp = t_exp;


[Nrows, Ncols] = size(t_exp);


for i=1:Ncols
     date_from_launch(i) = date(i) - launch_date;
     
     fprintf(filename,'%d %d %d \n', date(i), date_from_launch(i), texp(i));   
end

    
    fclose(filename);
return

