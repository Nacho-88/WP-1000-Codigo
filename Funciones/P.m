function P = P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la presión de la atmósfera a partir de la altitud (real/geométrica)
% Entrada: z (altitud geometrica)
% Salida: P (presion en Pascales)
% Si T_(n+1) ~= T_n -> P(z) = P_n*(1+[((z*(R_T/(R_T+z))-z_n)*(T_(n+1)-T_n))/((z_(n+1)-z_n)*T_n)])^((-g_0*(z_(n+1)-z_n))/(R*(T_(n+1)-T_n)))
% Si T_(n+1) = T_n -> P(z) = P_n*e^[(-g_0*(z*(R_T/(R_T+z))-z_n))/(R*T_n)]
% P_0 = 1 atm = 101325 Pa
% P_(n+1) = P(z_(n+1))
% Primero se introducen los valores de Z_n y T_n que son fijos, y luego se calculan todos los valores de P_n (lo ideal es guardarlos ya en vectores o matrices para seguir calculando)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constantes
R_T= 6371000;  % Radio de la Tierra medio (m)
g_0 = 9.81;      % Acelaración gravitatoria (m/s^2)
R  = 287.053;  % J/(kg·K)

% Vectores Altitud (geopotencial) y Temperatura (K)
z_n_geop = [0, 11000, 20000, 32000, 47000, 51000, 71000];
T_n = [288, 216.5, 216.5, 228.5, 270.5, 270.5, 214.5];

% Pasar de geopotencial a geometrica
z_n_geom = (R_T .* z_n_geop) ./ (R_T - z_n_geop);

% Calcular P_n capa a capa
P_n = generar_P_n();

% Localizar capa
n = n_capa(z);

% Evaluar P(z)

dz  = z_n_geom(n+1) - z_n_geom(n);       
dT  = T_n(n+1) - T_n(n);

if dT ~= 0
    P = P_n(n) * (1 + (((z * (R_T / (R_T + z)) - z_n_geom(n)) * dT) / (dz * T_n(n)))) ^ ((-g_0 * dz) / (R * dT));
else
    P = P_n(n) * exp((-g_0 * ((z * (R_T / (R_T + z)) - z_n_geom(n)))) / (R * T_n(n)));
end

end

