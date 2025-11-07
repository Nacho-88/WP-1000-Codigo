function T = T_atm(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la temperatura atmosférica estándar
% Entrada: z (altitud geometrica)
% Salida: T (temperatura en Kelvin)
% Formula Altura: z_geom = (R_t * z_geop) / (R_t - z_geop)
% Formula Temperatura: T_m(z) = T_n + (z - z_n_geom) * ((T_(n+1)-T_n) / (z_(n+1,geom)-z_(n,geom)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constante
R_T= 6371000; % Radio de la Tierra medio (m)

% Vectores Altitud (geopotencial) y Temperatura (K)
z_n_geop = [0, 11000, 20000, 32000, 47000, 51000, 71000];
T_n = [288, 216.5, 216.5, 228.5, 270.5, 270.5, 214.5];

% Pasar de geopotencial a geometrica
z_n_geom = (R_T .* z_n_geop) ./ (R_T - z_n_geop);

% Buscar la capa
n = find(z >= z_n_geom(1:end-1) & z < z_n_geom(2:end), 1, 'first');

if isempty(n)
    if z < z_n_geom(1)
        n = 1; % por debajo de 0
    else
        n = length(z_n_geom) - 1; % por encima de 71000
    end
end

% Calculo de la Tempearatura (K)

T = T_n(n) + (z - z_n_geom(n))*((T_n(n+1) - T_n(n)) / (z_n_geom(n+1) - z_n_geom(n)));

end

