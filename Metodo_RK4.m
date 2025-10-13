%Método de Runge-Kutta

function [times, y_aprox] = Metodo_RK4(dt, t0, tf, f, y0)

    times = t0 : dt : tf;

    Npts = length(times);

    % Número de EDOs del sistema
    n = length(y0);

    % Creamos el vector para guardar la solución aproximada
    y_aprox = zeros(n, Npts);
    y_aprox(:,1) = y0;

    % Resolvemos por el método de Runge-Kutta 4
    for k = 1 : Npts-1

        tk = times(k);
        yk = y_aprox(:,k);

        K1 = f(tk, yk);
        K2 = f(tk + dt/2, yk + (dt/2).*K1);
        K3 = f(tk + dt/2, yk + (dt/2).*K2);
        K4 = f(tk + dt, yk + dt.*K3);

        y_aprox(:,k+1) = yk + (dt/6).*(K1 + 2.*K2 + 2.*K3 + K4);
    
    end

end


