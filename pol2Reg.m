function [xmin, ymin, hess, it] = pol2Reg()

    xpoints = linspace(-20, 20);
    randos = rand([1,length(xpoints)]) - .5;
    function res = points(X)
        
       res = 3 + xpoints.*2 + xpoints.^2 + 60.*randos;
    end

    ypoints = points(xpoints);

    function y = curve(a, b, c, x)
       y = a + b.*x + c.*(x.^2); 
    end

    function r2 = func(a,b,c)
        r2 = 0;
        for i=1:length(xpoints)
            r2 = r2 + (ypoints(i) - a - b.*xpoints(i) - c.*(xpoints(i)^2)).^2;
        end
    end

    function fa = aGrad(a,b,c)
        fa = 0;
        for i = 1:length(xpoints)
            fa = fa + ypoints(i) - a - b.*xpoints(i) - c.*(xpoints(i)^2);
        end
        fa = -2*fa;
    end

    function fb = bGrad(a,b,c)
        fb = 0;
        for i = 1:length(xpoints)
            fb = fb + (ypoints(i) - a - b.*xpoints(i) - c.*(xpoints(i)^2)).*xpoints(i);
        end
        fb = -2*fb;
    end

    function fc = cGrad(a,b,c)
        fc = 0;
        for i = 1:length(xpoints)
            fc = fc + (ypoints(i) - a - b.*xpoints(i) - c.*(xpoints(i)^2)).*(xpoints(i).^2);
        end
        fc = -2*fc;
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
    [A,B,C] = meshgrid(-10:.1:10);
    Z = func(A,B,C);
    FA = aGrad(A,B,C);
    FB = bGrad(A,B,C);
    FC = cGrad(A,B,C);
    
    subplot(1,2,1)
    contour3(A(:,:,1), B(:,:,1), Z(:,:,1), 30)
    xlabel('A')
    ylabel('B')
    zlabel('Z')
    hold all
    n = 5;
    ai = 1.2;
    bi = 0;
    ci = 0;
    zi = func(ai, bi, ci);
    fa = 1;
    fb = 1;
    fc = 1;
    scatter3(ai, bi, zi, 'gs')
    lambdas = logspace(-8, 1, 50);
    it = 0;
    
%     for i = 1:n
    while abs([fa, fb, fc]) > .000001
        it  = it + 1;
        fa = aGrad(ai,bi, ci);
        fb = bGrad(ai,bi,ci);
        fc = cGrad(ai,bi,ci);
        Ap = ai - lambdas*fa;
        Bp = bi - lambdas*fb;
        Cp = ci - lambdas*fc;
        Zp = func(Ap, Bp, Cp);
        [~, ind] = min(Zp(1:46));
        lambda = lambdas(ind);
        af = ai - lambda*fa;
        bf = bi - lambda*fb;
        cf = ci - lambda*fc;
        zf = func(af, bf, cf);
        subplot(1,2,1)
        scatter3(af, bf, zf)
        plot3([ai,af], [bi,bf], [zi,zf])
        ai = af;
        bi = bf;
        ci = cf;
        zi = zf;
    end
%     hess = Hess(ai, bi);
    amin = ai
    bmin = bi
    cmin = ci
    subplot(1,2,2)
    scatter(xpoints, ypoints)
    xlabel('x')
    ylabel('y')
    hold all
    
    function res = plotter(a, b, c, x)
        res = a + b.*x + c.*(x.^2);
    end
    x = linspace(-20, 20);
    y = plotter(amin, bmin, cmin, x);
    plot(x, y)
    
    
    
end
