% to run the New Horizons pipeline, copy this file to a file named 'startup.m'
% and change the addpath commands below to match your file system.
% email zemcov@cfd.rit.edu with questions or comments.

%format short g
%format compact
dbstop if error
% more on

disp('[startup.m] reseeding random number generators')
rand('state',sum(100*clock))
randn('state',sum(100*clock))

%set(0,'DefaultAxesFontSize',12)
%set(0,'DefaultAxesFontName','Times')
% set(0,'DefaultLineMarkerSize',5)

% set(0,'DefaultFigurepaperpositionmode','auto')
% sets the default to print the figure as it appears in 
% the figure window

% ------------------------------------------
% add my matlab directories to path
% last directories listed are first to be searched!

addpath (genpath('/home/zemcov/NH/nh_analysis/matlab'))
addpath (genpath('/home/zemcov/matlab_lib'))

%%% end
