function nh_regress_dgl()

   load('../scratch/nh_make_results_nodgl.mat')

   bfac = (1. - 1.1 .* 0.75 .* sqrt(sin(pi.*[28.41,85.74,57.68,62.03]./180)))./0.35;
   
   ohm = ohm - 0.8;
   cobmean(3) = 2.*cobmean(3);
   
   ohmerr = sqrt(0.06.^2 + 0.027.^2).*ones(4,1);
   
   plot(ohm,cobmean,'o')
   errorbar(ohm,cobmean,coberr,coberr,ohmerr,ohmerr,'o')


   
   [a,b,alpha,p,chiopt,Cab,Calphap]=wtls_line(ohm,cobmean,ohmerr,coberr);
   
   brescal = a .* 100e-6 ./ (3e8 .* 1e-26 .* 1e6 .* 1e9)
   
   %disp(sprintf('Best-fitting
   
   disp(sprintf('p-value of fit: %f',1-chi2cdf(chiopt,2)))
   
   dbstop
   

% end