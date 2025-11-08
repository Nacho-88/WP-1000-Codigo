function T = T_globo(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la temperatura interna del globo a partir de la altitud (real/geométrica)
% Entrada: z (altitud geometrica)
% Salida: T (medida en Kelvin)
% z_star = 12000 m
% Depende de T_atm(z) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constantes
R_T= 6371000;   % Radio de la Tierra medio (m)
z_star = 12000; % Límite para cambio de formula (m)

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
T_k = generar_T_k(P_k);

% Zona Baja
if z <= z_star
    T = T_atm(z);
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
    Tk = T_k(k);

    P_z = P(z);

    T = Tk * (P_z / Pk)^(2/5);
end


end

