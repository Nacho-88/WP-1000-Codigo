function [r, a, b, eps_N_calc] = Met_Biseccion(f, a_0, b_0, param, eps_N)
% Hallar la raíz 'r' de la función 'f' en el intervalo cerrado [a,b] con
% una tolerancia eps y en N iteraciones. (Este método converge siempre que 'f' sea continua en el intervalo)
%
%-'f' es la función (anónima).
%
%-'a_0' y 'b_0' son los extremos del intervalo cerrado.
%
%-'param' es  un parámetro que indica si se introduce el valor de 'eps' (param=0) o 'N' (param=1).
%
%-'eps_N' es el error relativo máximo admitido en el valor de la raíz ('eps') o el
%número de iteraciones necesarios para la convergencia ('N') en función del valor de 'param'.
%
%
%-'r' es la raíz de la función en el intervalo.
%
%-'a' y 'b' son vectores (columna) tales que se almacenan los extremos inferior y
%superior del intervalo i-ésimo usado en las componentes a(i+1) y b(i+1). (Además a(1)=a_0 y b(1)=b_0)
%
%-'eps_N_calc' es el valor calculado de 'eps' o 'N' en función de cuál no
%se haya dado, es decir, si param=0 (eps_N=eps) entonces eps_N_calc=N, si
%param=1 (eps_N = N) entonces eps_N_calc=eps.

if param==0         %eps_N = eps
    eps = eps_N;
    N = ceil(log2((b_0-a_0)/eps));
    eps_N_calc = N;
elseif param==1     %eps_N = N
    N = eps_N;
    eps = (b_0-a_0)/2^N;
    eps_N_calc = eps;
else
    disp("Error en el valor de la variable 'param': Solo puede valer 1 o 0.")
    r = NaN;
    a=a_0;
    b=b_0;
    eps_N_calc=NaN;
    return
end

%Las variables 'eps' y 'N' tienen sus respectivos valores
%independientemente del valor de 'eps_N' y de 'param', por lo que usaremos
%estas variables auxiliares.

a = zeros(N, 1);
b = zeros(N, 1);

a(1,1) = a_0;
b(1,1) = b_0;

for i=1:N-1
        
    c = a(i,1) + (b(i,1)-a(i,1))/2;

    if sign(f(a(i,1)))~=sign(f(c))
        a(i+1,1) = a(i,1);
        b(i+1,1) = c;
    else
        a(i+1,1) = c;
        b(i+1,1) = b(i,1);
    end

end

r = a(end,1) + (b(end,1)-a(end,1))/2;

end

