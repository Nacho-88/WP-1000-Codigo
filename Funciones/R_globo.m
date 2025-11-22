function R_g = R_globo(z, m_He)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula el radio del globo a partir de la altitud (real/geométrica) suponiendo que mantiene forma de esfera
% Entrada: z (altitud geometrica), m_he (masa de Helio en kg)
% Salida: T (medida en Kelvin)
% Depende de T_atm(z) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Llamar constantes
load constantes R_prima

Tz = T_atm(z);
Pz = P(z);
R_g = ((3*m_He*R_prima)/(4*pi) * (Tz / Pz))^(1/3);

end

