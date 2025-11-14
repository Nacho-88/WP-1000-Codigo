function R_g = R_globo(z, m_He)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula el radio del globo a partir de la altitud (real/geométrica) suponiendo que mantiene forma de esfera
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg)
% Salida: T (medida en Kelvin)
% z_star = 12000 m
% Depende de T_atm(z) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Llamar constantes
load constantes R_prima

Tz = T_atm(z);
Pz = P(z);
R_g = ((3*m_He*R_prima)/(4*pi) * (Tz / Pz))^(1/3);

% % Determinar la capa (n*)
% n_star = n_capa(z);
% 
% % Definir los tramos del ascenso
% z_star_frontera = [z_star, z_n(n_star+1:end)];
% 
% % Zona Baja
%  if z <= z_star
%      Tz = T_globo(z);
%      Pz = P(z);
%      R_g = ((3*m_He*R_prima)/(4*pi) * (Tz / Pz))^(1/3);
%      return
%  end
% 
% % Zona Alta
% if z > z_star
%     % Determinar el tramo de la altitud
%     k = find(z >= z_star_frontera, 1, 'last');
%     if isempty(k)
%         k = 1;
%     end
% 
%     Pk = P_k(k);
%     Rk = R_k(k);
% 
%     P_z = P(z);
% 
%     R_g = Rk * (Pk / P_z)^(2/5);
% end

end

