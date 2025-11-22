function empuje = E(z, m_He)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la fuerza de empuje del globo que permite el ascenso a partir de la altitud (real/geométrica).
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg)
% Salida: empuje (medido en Newtons)
% Depende de R_globo(z,mHe) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Radio del globo
R_g = R_globo(z, m_He);

if R_g <= 0
   empuje = 0;
   warning('Radio de globo negativo, se toma como 0. (E)')
   return;
end


% Intervalos par para el método de Simpson
N = 20;

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

