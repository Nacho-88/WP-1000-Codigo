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


% Llamar constantes
load constantes R_T g_0 R z_n T_n P_n


% Localizar capa
n = n_capa(z);

% Evaluar P(z)

dz  = z_n(n+1) - z_n(n);       
dT  = T_n(n+1) - T_n(n);

if dT ~= 0
    P = P_n(n) * (1 + (((z * (R_T / (R_T + z)) - z_n(n)) * dT) / (dz * T_n(n)))) ^ ((-g_0 * dz) / (R * dT));
else
    P = P_n(n) * exp((-g_0 * ((z * (R_T / (R_T + z)) - z_n(n)))) / (R * T_n(n)));
end

end

