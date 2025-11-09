function integral = regla_de_simpson(dx, x, f)

    % fvals = f(x);

    fvals = zeros(size(x));
    for i=1:length(fvals)
        fvals(i) = f(x(i));
    end

    integral = (dx/3) * (fvals(1) + 4*sum(fvals(2 : 2 : end-1)) + 2*sum(fvals(3 : 2 : end-1)) + fvals(end));

end