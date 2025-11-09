function rho = rho_atm(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la densidad atmosférica del aire a partir de la altitud (real/geométrica)
% Entrada: z (altitud geometrica)
% Salida: rho (densidad en kg/m^3)
% rho_atm(z)=[(R*T_atm(z))/(P(z))]
% R ( J / (kg * K) ) = 287,053
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Llamar constantes
load constantes R

% Temperatura y Presión
Tz = T_atm(z);  % K
Pz = P(z);      % Pa

% Calcular densidad
rho = Pz./(R.*Tz);

end

