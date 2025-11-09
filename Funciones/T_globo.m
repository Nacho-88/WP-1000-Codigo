function T = T_globo(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la temperatura interna del globo a partir de la altitud (real/geométrica)
% Entrada: z (altitud geometrica)
% Salida: T (medida en Kelvin)
% z_star = 12000 m
% Depende de T_atm(z) y P(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Llamar constantes
load constantes z_n T_k z_star P_k

% Determinar la capa (n*)
n_star = n_capa(z);

% Definir los tramos del ascenso
z_star_frontera = [z_star, z_n(n_star+1:end)];

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

