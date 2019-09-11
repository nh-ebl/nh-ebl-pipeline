function [l,b]=RA_Dec_to_Gal(RA,Dec)
% RA_Dec_to_Gal: Convert RA and Dec to galactic coordinates
%
% [l,b]=RA_Dec_to_Gal(RA,Dec)
%
% ARGUMENTS
%  RA	The Right Ascension [degrees]
%  Dec	The Declination [degrees]
%
% RETURNS
%  l	Galactic longitude [degrees]
%  b	Galactic latitude [degrees]
%
% NOTES
%  Presuming J2000.
%
%  Uses the approximate formulae from "Allen's Astrophysical
%  Quantities", Aurthur N. Cox, Ed. 2000.
%  C.A. Murray, 1988, A&Ap 218, 325 supposedly has a more precise formula.
%
% See also: PRECESS

% AUTHOR: Eric Tittley
%
% HISTORY:
%  00 11 17 First version
%  07 11 05 Modified Comments
%
% COMPATIBILITY: Matlab, Octave
%
% LICENSE
%  Copyright Eric Tittley 2000
%  See the Gnu GPL license.
%  Essentially, free to use.  Free to modify.  But cannot be re-created in
%  whole or in part in anything sold or traded for something of worth.

Dec=Dec/180*pi;
RA=RA/180*pi;

Coef1=62.87/180*pi;
Coef2=282.86/180*pi;
Coef3=32.9319186/180*pi;
Coef4=2*pi+Coef3;

b=asin(sin(Dec)*cos(Coef1)-cos(Dec).*sin(RA-Coef2)*sin(Coef1));
l=0*b;

% Catch where cos(b)==0 and set the values by hand.  These correspond to the
% Galactic poles.
b_eq_2pi = find(cos(b)==0);
b(b_eq_2pi) = b(b_eq_2pi)*0+pi/2;
l(b_eq_2pi) = l(b_eq_2pi)*0;

% This is the break that separates into hemispheres the domains for which we
% have to get the signs right.  We also catch cos(b)==0, which would give a
% divide by zero error.
Dec_le_break = find( (Dec <= (Coef1-pi/2)*sin(RA-Coef2)) & cos(b)~=0 );
Dec_gt_break = find( (Dec  > (Coef1-pi/2)*sin(RA-Coef2)) & cos(b)~=0 );

l(Dec_gt_break)=Coef3 + ...
 acos( cos(Dec(Dec_gt_break)).*cos(RA(Dec_gt_break)-Coef2)./cos(b(Dec_gt_break)));

l(Dec_le_break)=Coef4 - ...
 acos( cos(Dec(Dec_le_break)).*cos(RA(Dec_le_break)-Coef2)./cos(b(Dec_le_break)));

l_ge_2pi = find(l>=2*pi);
l(l_ge_2pi)=l(l_ge_2pi)-2*pi;

b=b*180/pi;
l=l*180/pi;