function R_g = R_globo(z, m_he)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula el radio del globo a partir de la altitud (real/geométrica) suponiendo que mantiene forma de esfera
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg)
% Salida: T (medida en Kelvin)
% z_star = 12000 m
% Depende de T_atm(z) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constantes
R_T = 6371000;     % Radio medio Tierra (m)
z_star = 12000;    % Altura de transición (m)

R = 287.053;          % J/(kg*K)
m_a   = 0.0289644;        % kg/mol
m_He_molar = 0.0040026;   % kg/mol

% R_prima = R'(apuntes)
R_prima = (m_a / m_He_molar) * R; % J/(kg*K)

% Vector Altitud (geopotencial)
z_n_geop = [0, 11000, 20000, 32000, 47000, 51000, 71000];

% Pasar de geopotencial a geometrica
z_n_geom = (R_T .* z_n_geop) ./ (R_T - z_n_geop);

% Determinar la capa (n*)
n_star = n_capa(z);

% Definir los tramos del ascenso
z_star_frontera = [z_star, z_n_geom(n_star+1:end)];

% Llamar P_k y T_k
P_k = generar_P_k();          
R_k = generar_R_k(m_he);

% Zona Baja
 if z <= z_star
     Tz = T_globo(z);
     Pz = P(z);
     R_g = ((3*m_he*R_prima)/(4*pi) * (Tz / Pz))^(1/3);
     return
 end

% Zona Alta
if z > z_star
    % Determinar el tramo de la altitud
    k = find(z >= z_star_frontera, 1, 'last');
    if isempty(k)
        k = 1;
    end

    Pk = P_k(k);
    Rk = R_k(k);

    P_z = P(z);

    R_g = Rk * (Pk / P_z)^(2/5);
end

end

