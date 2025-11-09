function F = F_roz_globo(z, dz_dt, m_He)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la la fuerza de rozamiento del globo a partir de la altitud (real) y su velocidad
% Entrada: z (altitud geometrica), dz_dt (velocidad ascenso en m/s)
% Salida: Fuerza: F (medida en Newtons)
% Constante: C_D_semiesfera
% Depende del Radio del globo (R_globo) y Densidad atmosfera (rho_Atm), Presión atmósfera (P) y Integración (Método Simpson 1/3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Llamar constantes
load constantes C_D_globo

% Radio del globo
R_z = R_globo(z, m_He);

if R_z <= 0
   F = 0;
   warning('Radio de globo negativo, se toma como 0. (F_roz_globo)')
   return;
end

% Intervalos par para el método de Simpson
N = 200;

if mod(N,2) ~= 0
   N = N + 1;
end

% Integración
a = 0;
b = R_z;

x = linspace(a, b, N+1);  % vector de nodos
dx = (b - a) / N;         % paso entre nodos
f = @(zz) rho_atm(z + zz) .* zz;

I_h = regla_de_simpson(dx, x, f);

% Fuerza de rozamiento
F = -pi*(dz_dt.^2)*C_D_globo*I_h;

end

