% MATLAB script to reduce ra/dec data by n degrees
% This will seperate the data into feilds, such that catalog data retrival
% is easier. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pgi8114
% 2/15/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ra_dec_reduce(phase)

% load matrix with RA_array, DEC_array, date_array, file_array.
load(sprintf('../../data/NH/matrix/%s_ra_dec.mat',phase));
import matlab.io.*

DEC = DEC_array;
RA = RA_array;

% define arrays
RA_reduced=[]; %in degrees
DEC_reduced=[]; %in degress
%RA_sexagesmal = [];
%sexagesmal = [];
RAhrs_reduced = []; %RA in hours
start_index=[];
end_index=[];


start_index = [start_index 1];
RA_reduced = [RA_reduced RA(start_index)];
DEC_reduced = [DEC_reduced DEC(start_index)];

arcminutes = 30; %maximun value given by catalog 
degrees = arcminutes/60 - .07; % .07 is the maximun distance in degrees away
                                %from the center that a pix can be.
n = size(RA,2);
j=1;

for curr_index = 2:n
	if (RA(curr_index) > RA_reduced(j) + degrees) || (RA(curr_index) < RA_reduced(j) - degrees) || (DEC(curr_index) > DEC_reduced(j) + degrees) || (DEC(curr_index) < DEC_reduced(j) - degrees)
		end_index = [end_index curr_index];
		j = j+1;
		if curr_index < n
			start_index = [start_index curr_index+1];
			RA_reduced = [RA_reduced RA(curr_index+1)];
			DEC_reduced = [DEC_reduced DEC(curr_index+1)];
		end
	end
end

if size(start_index,2) > size(end_index,2)
	end_index = [end_index n];
end

% Convert RA to hours
[m,n] = size(RA_reduced);
for i = 1:n
	RAcurr = RA_reduced(i)/15;
	RAhrs_reduced = [RAhrs_reduced RAcurr];
end

%%convert RA and Dec to sexigesimal format.
%[m n] = size(RA_reduced);
%%To use, fix DEC conversion for negatives
%for i = 1:n
%	HRS = RA_reduced(i)/15;
%	RAh = floor(HRS);
%	MIN = (HRS - RAh)*60;
%	RAm = floor(MIN);
%	RAs = (MIN - RAm)*60;
%	
%	%RAcurr = sprintf('%02.0f:%02.0f:%02.1f \n');
%	%RA_sexagesmal = [RA_sexagesmal RAcurr];
%	
%	DEG = DEC_reduced(i);
%	DECh = floor(DEG);
%	MIN = (DEG - DECh)*60;
%	DECm = floor(MIN);
%	DECs = (MIN - DECm)*60;
%
%	curr = sprintf('%02.0f:%02.0f:%02.1f  %02.0f:%02.0f:%02.1f \n', RAh, RAm, RAs, DECh, DECm, DECs);
%	sexagesmal = [sexagesmal curr];
%end

%save
save(sprintf('../../data/NH/matrix/%s_ra_dec_reduced.mat', phase), 'RA_reduced','DEC_reduced', 'RAhrs_reduced', 'start_index', 'end_index');

return
