%-------------------------------------------------------------------------
%
% Planpos: Ephemerides of the major planets
%
% Reference:
% Montenbruck O., Pfleger T.; Astronomy on the Personal Computer; Springer 
% Verlag, Heidelberg; 4th edition (2000).
%
% Last modified:   2015/08/12   M. Mahooti
% Adapted to be used as a function: 2019/04/22 T. Symons
%-------------------------------------------------------------------------
function planet_pos(Mode,year,month,day,Hour)
Ast_Const

End = false;

% Command loop
while(1)
   
    % Epoch date
    fprintf(' Date (yyyy mm dd hh.hhh) ... \n');
    year = input(' yyyy ');
    month = input(' mm ');
    day = input(' dd ');
    Hour = input(' hh.hhh ');
    fprintf('\n');
    MJD = Mjday(year,month,day)+Hour/24;
    T   = (MJD-MJD_J2000)/36525;    
    % Header
    [year,mon,day,hour,minute,sec] = invjday(MJD+2400000.5);
    
    % Ecliptic coordinates of the Sun, equinox of date
    R_Sun = SunPos(T);    
    % Loop over planets
    if ( (-1.1<T) && (T<1) )
        last = 10;
    else
        last = 9;
    end    
    for iPlanet = 1:last
        % Heliocentric ecliptic coordinates of the planet; equinox of date
        switch iPlanet
            case 1
                Planet = 'Sun';
            case 2
                Planet = 'Mercury';
            case 3
                Planet = 'Venus';
            case 4
                Planet = 'Earth';
            case 5
                Planet = 'Mars';
            case 6
                Planet = 'Jupiter';
            case 7
                Planet = 'Saturn';
            case 8
                Planet = 'Uranus';
            case 9
                Planet = 'Neptune';
            case 10
                Planet = 'Pluto';
        end        
        r_helioc = PertPosition(Planet, T);        
        % Geocentric ecliptic coordinates (equinox of date)
        r_geoc = r_helioc + R_Sun;        
        % First-order light-time/aberration correction
        dist = norm(r_geoc);
        fac = dist/c_light;        
        if (strcmp(Mode,'a'))
            r_geoc = r_geoc - fac*(KepVelocity(Planet,T)-KepVelocity('Earth',T));
        else
            r_geoc = r_geoc - fac*KepVelocity(Planet,T);
        end        
        % Precession and equatorial coordinates
        switch Mode
            case 'a'
                r_equ    = NutMatrix(T) * Ecl2EquMatrix(T) * r_geoc;
            case 'j'
                P        = PrecMatrix_Ecl(T,T_J2000);
                r_helioc = P * r_helioc;
                r_geoc   = P * r_geoc;
                r_equ    = Ecl2EquMatrix(T_J2000) * r_geoc;
            case 'b'
                P        = PrecMatrix_Ecl(T,T_B1950);
                r_helioc = P * r_helioc;
                r_geoc   = P * r_geoc;
                r_equ    = Ecl2EquMatrix(T_B1950) * r_geoc;
        end        
        % Output
        fprintf(' %8s', Planet);
        [r_helioc_phi, r_helioc_theta, r_helioc_r] = CalcPolarAngles(r_helioc);
        [D, M, S] = DMS(Deg*r_helioc_phi);
        fprintf('%4d%3.2d%5.1f ', D, M, S);
        [D, M, S] = DMS(Deg*r_helioc_theta);
        fprintf('%4d%3.2d%5.1f ', D, M, S);
        if (iPlanet<=4)
            fprintf('%11.6f', r_helioc_r);
        else
            fprintf('%10.5f ', r_helioc_r);
        end
        [r_equ_phi, r_equ_theta, r_equ_r] = CalcPolarAngles(r_equ);
        [D, M, S] = DMS(Deg*r_equ_phi/15);
        fprintf('%4d%3.2d%6.2f ', D, M, S);
        [D, M, S] = DMS(Deg*r_equ_theta);
        fprintf('%4d%3.2d%5.1f ', D, M, S);
        if (iPlanet<=4)
            fprintf('%11.6f', dist);
        else
            fprintf('%10.5f ', dist);
        end
        fprintf('\n');
    end    
    % Trailer
    fprintf('\n l,b,r:   heliocentric ecliptic (geometric) \n');
    fprintf(' RA,Dec:  geocentric equatorial ');
    if (strcmp(Mode,'a'))
        fprintf('(apparent)\n');
    else
        fprintf('(astrometric)\n');
    end
    fprintf(' delta:   geocentric distance   (geometric)\n\n');    
    if (~End)
        break
    end
end
end

