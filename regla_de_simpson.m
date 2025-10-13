function integral = regla_de_simpson(dx, x, f)

    fvals = f(x);

    integral = (dx/3) * (fvals(1) + 4*sum(fvals(2 : 2 : end-1)) + 2*sum(fvals(3 : 2 : end-1)) + fvals(end));

end