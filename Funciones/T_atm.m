function T = T_atm(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la temperatura atmosférica estándar
% Entrada: z (altitud geometrica)
% Salida: T (temperatura en Kelvin)
% Formula Altura: z_geom = (R_t * z_geop) / (R_t - z_geop)
% Formula Temperatura: T_m(z) = T_n + (z - z_n_geom) * ((T_(n+1)-T_n) / (z_(n+1,geom)-z_(n,geom)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Llamar constantes
load constantes z_n T_n R_T

% Buscar la capa
n = n_capa(z);

% Calculo de la Tempearatura (K)
T = T_n(n) + ((z*(R_T/(R_T+z)))-z_n(n))*((T_n(n+1)-T_n(n))/(z_n(n+1)-z_n(n)));

end

