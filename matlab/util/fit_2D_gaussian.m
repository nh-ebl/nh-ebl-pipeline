function [fitresult, gof] = fit_2D_gaussian(x, y, z, w)
%fit_2D_Gaussian(X,Y,Z)
%  Fit 2D Gaussian to data.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xx
%      Y Input : yy
%      Z output: zz
%      Weights : w
%
%  Output:
%      fitresult : an sfit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, SFIT.

%% Initialization.

%comment
% Convert all inputs to column vectors.
z(~isfinite(z)) = 0;
xx = x(:);
yy = y(:);
zz = z(:);

%% Fit: '2D Gaussian'.
ft = fittype( 'A*exp(-0.5*((x-xo)^2/(fwhm_x/2.35482)^2+(y-yo)^2/(fwhm_y/2.35482)^2))+C', 'indep', {'x', 'y'}, 'depend', 'z' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf 0 0 0 0];
opts.StartPoint = [1 0.1 20 20 50 50];
opts.Upper = [Inf Inf Inf Inf 1024 1024];
opts.Weights = zeros(1,0);
[fitresult, gof] = fit( [xx, yy], zz, ft, opts );


% % Plot fit with data.
% figure( 'Name', '2D Gaussian Fit' );
% h = plot( fitresult, [xx, yy], zz );
% legend( h, '2D Gaussian Fit', 'zz vs. xx, yy', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'xx' );
% ylabel( 'yy' );
% zlabel( 'zz' );
%grid on


