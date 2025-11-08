function empuje = E(z, m_he)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la fuerza de empuje del globo que permite el ascenso a partir de la altitud (real/geométrica).
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg)
% Salida: empuje (medido en Newtons)
% z_star = 12000 m
% Depende de R_globo(z,mHe) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Radio del globo
R_g = R_globo(z, m_he);

if R_g <= 0
   empuje = 0;
   return;
end

% Intervalos par para el método de Simpson
N = 200;

if mod(N,2) ~= 0
  N = N + 1;
end

% Integración
a = -R_g;
b = R_g;

x = linspace(a, b, N+1);  % vector de nodos
dx = (b - a) / N;         % paso entre nodos
f = @(zz) P(z + zz) .* zz;

I_h = regla_de_simpson(dx, x, f);

% Empuje
empuje = -2*pi * I_h;

end

