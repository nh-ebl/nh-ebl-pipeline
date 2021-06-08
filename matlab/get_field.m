function [ field ] = get_field(ra,dec,paths)
% Get the field number for a fits file.
%   Returns the number of the corresponding catalog field based on the ra 
%   and dec of the input file. The values for each field are based on the
%   output values of ra_dec_reduce.m and the catalog fields were manully
%   retrieved and converted to matrices using catalog_data.m
%
% PARAMETERS:
%       ra,dec = ra,dec of the field to be enumerated (in degrees).
%
% OUTPUT:
%        field = the number of the correspoding catalog field based on the
%                ra and dec of the input file. 
%
% author: Poppy Immel
% date: May 3, 2016
% email: pgi8114@rit.edu
% modified: MZ, Jun 16 2016

  load(sprintf('%scataloginfo.mat',paths.catdir)); 
  [~,n] = size(field_number);

  arcminutes = 15; %maximum value given by catalog 
  degrees = arcminutes/60 - .07; 

  for i = 1:n
    if ra < field_RA(i) + degrees  && ra > field_RA(i) - degrees && dec < field_DEC(i) + degrees  && dec > field_DEC(i) - degrees
        field = field_number(i);
    end
  end
  
  if isempty(field)
    warn('There is no corresponding field.');
    field = NaN;
  end
    
end

