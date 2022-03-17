function [a, b, siga, sigb, chi2, q] = fitexy(xx, yy, sx, sy)
    %%%
    % Straight-line fit to input data x[1..ndat] and y[1..ndat] with errors in both x and y, the respective standard deviations being the input quantities sx[1..ndat] and sy[1..ndat].
    % Output quantities are a and b such that y = a + bx minimizes ?2, whose value is returned
    % as chi2. The ?2 probability is returned as q, a small value indicating a poor fit (sometimes
    % indicating underestimated errors). Standard errors on a and b are returned as siga and sigb.
    % These are not meaningful if either (i) the fit is poor, or (ii) b is so large that the data are
    % consistent with a vertical (infinite b) line. If siga and sigb are returned as BIG, then the data
    % are consistent with all values of b.
    %%%
    ndata = length(xx); %get size of data
    varx = var(xx); %get variances
    vary = var(yy);
    scale = sqrt(varx./vary); %get scaler
    yy = yy.*scale; %scale yy
    sy = sy.*scale; %scale sy
    ww = sqrt(sx.^2+sy.^2); %get weights

    [~, b, ~, ~, ~, ~] = fit_linearLeastSquares(xx, yy, ww); %linear least squares fit, get initial fit for b

    %get reference angles, and make b an angle
    offs = 0;
    ch = zeros(6,1);
    ang = zeros(6,1);
    ang(2) = atan(b);
    ang(5) = ang(2);
    ang(6) = pi/2;
    for j = 4:1:6
        [ch(j), ~] = chixy(xx, yy, sx, sy, offs, ang(j));
    end
    %Bracket the chiSq minimum
    [ang(1), ang(2), ang(3), ch(1), ch(2), ch(3)] = mnbrak(ang(1), ang(2), ang(3), xx, yy, sx, sy, offs);
    % then locate it with brent
    [~, b] = brent(ang(1), ang(2), ang(3), xx, yy, sx, sy, offs, 1e-3);
    [chi2, aa] = chixy(xx, yy, sx, sy, offs, b);
    a = aa;
    %compute chiSq probability
    q = gammainc(0.5*chi2,0.5*(ndata-2),'upper'); %gammq returns the incomplete gamma function Q(a,x) = 1-P(a,x) which is equivalent to the upper incomplete gamma function gammainc(x,a,'upper') - note x and a switched
    %save the inverse sum of weights at the minimum
    r2 = 1/sum(ww);
    bmx = 1e30;
    bmn = 1e30;
    offs = chi2+1;
    % now find standard errors for b as points where delta-chiSq = 1
    for j = 1:1:6 %go through saved values to bracket the desired roots. note periodicity in slope angles
        if ch(j) > offs
            d1 = abs(ang(j) - b);
            while d1 >= pi
                d1 = d1 - pi;
            end
            d2 = pi-d1;
            if ang(j) < b
                swap = d1;
                d1 = d2;
                d2 = swap;
            end
            if d1 < bmx
                bmx = d1;
            end
            if d2 < bmn
                bmn = d2;
            end
        end
    end
    if bmx < 1e30 % call zbrent to find the roots
        bmx = zbrent(b, b+bmx, 1e-3, xx, yy, sx, sy, offs) - b;
        amx = aa-a;
        bmn = zbrent(b, b-bmn, 1e-3, xx, yy, sx, sy, offs) - b;
        amn = aa-a;
        sigb = sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*(cos(b)^2));
        siga = sqrt(0.5*(amx*amx+amn*amn)+r2)/scale; %error has an additional piece r2
    else
        sigb = 1e30;
        siga = 1e30;
    end
    a = a/scale; %unscale answers
    b = tan(b)/scale;
end

function [a, b, siga, sigb, chi2, q] = fit_linearLeastSquares(xx, yy, sig)
    %%%
    % Adapted from Numerical Recipes in C book
    % Given a set of data points xx[1..ndata],yy[1..ndata] with individual standard deviations
    % sig[1..ndata], fit them to a straight line y = a + bx by minimizing ?2. Returned are
    % a,b and their respective probable uncertainties siga and sigb, the chi-square chi2, and the
    % goodness-of-fit probability q (that the fit would have ?2 this large or larger). If sig is not provided on
    % input, then the standard deviations are assumed to be unavailable: q is returned as 1.0 and
    % the normalization of chi2 is to unit standard deviation on all points.
    %%%
    ndata = length(xx);
    %sig is optional
    %accumulate sums
    if exist('sig','var')
        wt = 1./(sig.^2); %transform into weights
        ss = sum(wt); %sum of weights
        sx = sum(xx.*wt); %sum of weighted xx and yy
        sy = sum(yy.*wt);
        sxoss = sx/ss; %dunno
        t = (xx-sxoss)./sig;
        st2 = sum(t.^2);
        b = sum(t.*yy./sig);
    else
        sx = sum(xx); %no weighting
        sy = sum(yy);
        ss = length(xx); %size of xx is ss here w/o weights
        sxoss = sx/ss; %dunno
        t = x-sxoss;
        st2 = sum(t.^2);
        b = sum(t.*yy);
    end
    %solve for a, b, siga, sigb
    b = b/st2; 
    a = (sy-sx*b)/ss;
    siga = sqrt((1+sx^2/(ss*st2))/ss);
    sigb = sqrt(1/st2);
    if exist('sig','var')
        chi2 = sum(((yy-a-b.*xx)./sig).^2);
        if ndata > 2
            q = gammainc(0.5*chi2,0.5*(ndata-2),'upper'); %gammq returns the incomplete gamma function Q(a,x) = 1-P(a,x) which is equivalent to the upper incomplete gamma function gammainc(x,a,'upper') - note x and a switched
        else
            q = 1;
        end
    else
        chi2 = sum((yy-a-b.*xx).^2);
        sigdat = sqrt(chi2/(ndata-2));
        siga = siga*sigdat;
        sigb = sigb*sigdat;
        q = 1;
    end
    
end

function [answr, aa] = chixy(xx, yy, sx, sy, offs, bang)
    %%%
    % Adapted from Numerical Recipes in C book
    % Captive function of fitexy, returns the value of (?2 ? offs) for the slope b=tan(bang).
    % Scaled data and offs are NOT communicated via the global variables.
    %%%
    b = tan(bang);
    ww = (b.*sx).^2 + (sy).^2;
    k = ww < 1E-30;
    ww( k ) = 1E30; %seems to cap a 1/ww thing
    ww( ~k ) = 1./ww( ~k ); %flip the rest
    sumw = sum(ww);
    avex = sum(ww.*xx)/sumw;
    avey = sum(ww.*yy)/sumw;
    aa = avey - b.*avex;
    answr = sum(ww.*(yy-aa-b.*xx).^2) -offs;
end

function [ax, bx, cx, fa, fb, fc] = mnbrak(ax, bx, ~, xx, yy, sx, sy, offs) %~ is cx, really passed b/c of C pointer shennanigans
    %%%
    % Adapted from Numerical Recipes in C book
    % Given a function chixy, and given distinct initial points ax and bx, this routine searches in
    % the downhill direction (defined by the function as evaluated at the initial points) and returns
    % new points ax, bx, cx that bracket a minimum of the function. Also returned are the function
    % values at the three points, fa, fb, and fc.
    %%%
    gold = ( 1 + sqrt(5) ) / 2;
    glimit = 100;
    tiny = 1e-20;
    function [a, b, c] = shifty(~, b, c, d)
        a = b;
        b = c;
        c = d;
    end
    function return_a = sign_dub(a,b)
        if b >= 0
            return_a = abs(a);
        else
            return_a = -abs(a);
        end
    end

    fa = chixy(xx, yy, sx, sy, offs, ax);
    fb = chixy(xx, yy, sx, sy, offs, bx);
    if fb > fa
        % switch roles of a and b so that we can go downhill int he direction from a to b
        tmp = ax;
        ax = bx;
        bx = tmp;
        tmp = fa;
        fa = fb;
        fb = tmp;
    end
    %first guess for c
    cx = bx+gold*(bx-ax);
    fc = chixy(xx, yy, sx, sy, offs, cx);
    % keep returning here until we bracket
    %compute u by parabolic exptrapolation from a, b, c,
    % tiny is used to prevent any possible division by zero
    while fb > fc
        r = (bx-ax)*(fb-fc);
        q = (bx-cx)*(fb-fa);
        u = bx - ((bx-cx)*q - (bx-ax)*r)/(2*sign_dub(max(abs(q-r), tiny), q-r));
        ulim = bx + glimit*(cx-bx);
        %??we won't go farther than this. test various possibilities:
        if (bx-u)*(u-cx) > 0 %parabolic u is between b and c
            fu = chixy(xx, yy, sx, sy, offs, u);
            if fu < fc %got a minimum between b and c
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return
            elseif fu > fb %got a minimum between a and u
                cx = u;
                fc = fu;
                return
            else %parabolic fit was no use. use default magnifcation
                u = cx + gold*(cx-bx);
                fu = chixy(xx, yy, sx, sy, offs, u);
            end
        elseif (cx-u)*(u-ulim) > 0 %parabolic fit is between c and its allowed limit
            fu = chixy(xx, yy, sx, sy, offs, u);
            if fu < fc
                [bx, cx, u] = shifty(bx, cx, u, cx+gold*(cx-bx));
                [fb, fc, fu] = shifty(fb, fc, fu, chixy(xx, yy, sx, sy, offs, u));
            end
        elseif (u-ulim)*(ulim-cx) >= 0 %limit parabolic u to maximum allowed value
            u = ulim;
            fu = chixy(xx, yy, sx, sy, offs, u);
        else %reject parabolic u, use default magnetication
            u = cx + gold*(cx-bx);
            fu = chixy(xx, yy, sx, sy, offs, u);            
        end
        %eliminate oldest point and continue
        [ax, bx, cx] = shifty(ax,bx,cx,u);
        [fa, fb, fc] = shifty(fa,fb,fc,fu);
    end
end

function [fx, xmin] = brent(ax, bx, cx, xx, yy, sx, sy, offs, tol)
    %%%
    % Adapted from Numerical Recipes in C book
    % Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
    % between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
    % the minimum to a fractional precision of about tol using Brent’s method. The abscissa of
    % the minimum is returned as xmin, and the minimum function value is returned as brent, the
    % returned function value.
    %%%
    itmax = 100;
    cgold = 1-(( 1 + sqrt(5) ) / 2 -1); %some golden ratio thing idk
    zeps = 1e-10;
    % Here ITMAX is the maximum allowed number of iterations; CGOLD is the golden ratio; ZEPS is
    % a small number that protects against trying to achieve fractional accuracy for a minimum that
    % happens to be exactly zero.
    function [a, b, c] = shifty(~, b, c, d)
        a = b;
        b = c;
        c = d;
    end
    function return_a = sign_dub(a,b)
        if b >= 0
            return_a = abs(a);
        else
            return_a = -abs(a);
        end
    end
    
    e = 0; %this will be the distance moved on the step before last
    if ax < cx % a and b must be in ascending order, the input abscissas need not be
        a = ax;
        b = cx;
    else
        a = cx;
        b = ax;
    end
    x = bx; %init
    w = bx;
    v = bx;
    fw = chixy(xx, yy, sx, sy, offs, x);
    fv = chixy(xx, yy, sx, sy, offs, x);
    fx = chixy(xx, yy, sx, sy, offs, x);
    for iter = 1:itmax
        xm = 0.5*(a+b);
        tol1 = tol*abs(x)+zeps;
        tol2 = 2*tol1;
        if abs(x-xm) <= (tol2-0.5*(b-a)) %test for done here
            xmin = x;
            return 
        end
        if abs(e) > tol1 %construct a trial parabolic fit
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q - (x-w)*r;
            q = w*(q-r);
            if q > 0
                p = -p;
            end
            q = abs(q);
            etemp = e;
            e = d;
            if (abs(p) >= abs(0.5*(q*etemp))) || (p <= q*(a-x)) || (p >= q*(b-x))
                if x >= xm
                    e = a-x;
                else
                    e = b-x;
                end
                d = cgold*e;
                % The above conditions determine the acceptability of the parabolic fit. Here we
                % take the golden section step into the larger of the two segments.
            else
                %take the parabolic step
                d = p/q;
                u = x+d;
                if (u-a < tol2) || (b-u < tol2)
                    d = sign_dub(tol1, xm-x);
                end
            end
        else
            if x >= xm
                e = a-x;
            else
                e = b-x;
            end
            d = cgold*e;
        end
        if abs(d) >= tol1
            u = x+d;
        else
            u = x+sign_dub(tol1,d);
        end
        fu = chixy(xx, yy, sx, sy, offs, u);
        %this is one function evaluation per iteration
        if fu <= fx %now decide hwat to do with our function evaluation
            if u >= x
                a = x;
            else
                b = x;
            end
            [v, w, x] = shifty(v, w, x, u); %housekeeping follows
            [fv, fw, fx] = shifty(fv, fw, fx, fu);
        else
            if u < x
                a = u;
            else
                b = u;
            end
            if (fu <= fw) || (w == x)
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            elseif (fu <= fv) || (v == x) || (v == w)
                v = u;
                fv = fu;
            end
        end
        %done with housekeeping. back for another iteration
    end
    disp('Error: too many iterations in brent');
    xmin = x;
end

function b = zbrent(x1, x2, tol, xx, yy, sx, sy, offs)
    %%%
    % Adapted from Numerical Recipes in C book
    % Using Brent’s method, find the root of a function chixy known to lie between x1 and x2. The
    % root, returned as zbrent, will be refined until its accuracy is tol.
    %%%
    itmax = 100; %Maximum allowed number of iterations.
%     eps = 3e-8; %machine floating-point precision, f64 is higher and eps works in matlab already
    function return_a = sign_dub(a,b)
        if b >= 0
            return_a = abs(a);
        else
            return_a = -abs(a);
        end
    end

    a = x1;
    b = x2;
    c = x2;
    fa = chixy(xx, yy, sx, sy, offs, a);
    fb = chixy(xx, yy, sx, sy, offs, b);
    if ((fa > 0) && (fb > 0)) || ((fa < 0) && (fb < 0))
        error('Error in zbrent: root must be bracketed in zbrent');
    end
    fc = fb;
    for iter = 1:1:itmax
        disp(iter)
        if ((fb > 0) && (fc > 0)) || ((fa < 0) && (fb < 0))
            % rename a, b, c and adjust bounting interval d
            c = a;
            fc = fa;
            d = b - a;
            e = b - a;
        end
        if abs(fc) < abs(fb)
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        end
        tol1 = 2*eps*abs(b)+0.5*tol; %convergence check
        xm = 0.5*(c-b);
        if (abs(xm) <= tol1) || (fb == 0)
            return; %reuturn b
        end
        if (abs(e) >= tol1) && (abs(fa) > abs(fb))
            %attempt inverse quadratic interpolation
            s = fb/fa;
            if a == c
                p = 2*xm*s;
                q = 1-s;
            else
                q = fa/fc;
                r = fb/fc;
                p = s*(2*xm*q*(q-r)-(b-a)*(r-1));
                q = (q-1)*(r-1)*(s-1);
            end
            %check whether in bounds
            if p > 0
                q = -q;
            end
            p = abs(p);
            min1 = 3*xm*q-abs(tol1*q);
            min2 = abs(e*q);
            if min1 < min2
                test_min = min1;
            else
                test_min = min2;
            end
            if 2*p < test_min
                %accept interpolation
                e = d;
                d = p/q;
            else
                %interpolation failed, use bisection
                d = xm;
                e = d;
            end
        else
            %bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
        end
        a = b; %move last best guess to a
        fa = fb;
        if abs(d) > tol1 %evaluate new trial root
            b = b + d;
        else
            b = b + sign_dub(tol1,xm);
        end
        fb = chixy(xx, yy, sx, sy, offs, b);
    end
    disp('Maximum number of iterations exceeded in zbrent');
    b = 0;
end


