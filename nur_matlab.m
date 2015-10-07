function [xmin, ymin, hess, it] = linReg()

    xpoints = [1, 2, 3, 4, 5, 3, 7, 8, 9, 7];
    ypoints = [1, 4, 3, 2, 5, 6, 2, 8, 9, 10];

    function y = line(a, b, x)
       y = a + b.*x; 
    end
    function r2 = func(a, b)
        r2 = 0;
        for i=1:length(xpoints)
            r2 = r2 + (ypoints(i) - a - b.*xpoints(i)).^2;
        end
    end

    function fa = aGrad(a, b)
        fa = 0;
        for i = 1:length(xpoints)
            fa = fa + ypoints(i) - a - b.*xpoints(i);
        end
        fa = -2*fa;
    end

    function fb = bGrad(a, b)
        fb = 0;
        for i = 1:length(xpoints)
            fb = fb + (ypoints(i) - a - b.*xpoints(i)).*xpoints(i);
        end
        fb = -2*fb;
    end

    function fbb = Hes1(a,b)
        fbb = 0;
        for i = 1:length(xpoints)
            fbb = fbb + xpoints(i)^2;
        end
    end

    function faa = Hes2(a,b)
        faa = 2*length(xpoints);
    end

    function fab = Hes3(a,b)
        fab = -2*sum(xpoints);
    end

    function hess = Hess(a,b)
        fbb = Hes1(a,b);
        faa = Hes2(a,b);
        fab = Hes3(a,b);
        hess = fbb*faa - fab^2;
    end
    [A,B] = meshgrid(-10:.1:10, -10:.1:10);
    Z = func(A,B);
    FA = aGrad(A,B);
    FB = bGrad(A,B);
    
    subplot(1,2,1)
    contour3(A, B, Z, 30)
    xlabel('A')
    ylabel('B')
    zlabel('Z')
    hold all
    n = 5;
    ai = 1.2;
    bi = 0;
    zi = func(ai, bi);
    fa = 1;
    fb = 1;
    scatter3(ai, bi, zi, 'gs')
    lambdas = logspace(-8, 1, 50);
    it = 0;
    
%     for i = 1:n
    while abs([fa, fb]) > .00001
        it  = it + 1;
        fa = aGrad(ai, bi);
        fb = bGrad(ai, bi);
        Xp = ai - lambdas*fa;
        Yp = bi - lambdas*fb;
        Zp = func(Xp, Yp);
        [~, ind] = min(Zp(1:46));
        lambda = lambdas(ind);
        af = ai - lambda*fa;
        bf = bi - lambda*fb;
        zf = func(af, bf);
        subplot(1,2,1)
        scatter3(af, bf, zf)
        plot3([ai,af], [bi,bf], [zi,zf])
        ai = af;
        bi = bf;
        zi = zf;
    end
    hess = Hess(ai, bi);
    xmin = ai;
    ymin = bi;
    subplot(1,2,2)
    scatter(xpoints, ypoints)
    xlabel('x')
    ylabel('y')
    hold all
    plot([0, 10], [line(xmin, ymin, 0), line(xmin, ymin, 10)])
end
