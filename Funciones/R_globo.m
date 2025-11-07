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

% R_eff = R'(apuntes)
R_eff = (m_a / m_He_molar) * R; % J/(kg*K)

% Vector Altitud (geopotencial)
z_n_geop = [0, 11000, 20000, 32000, 47000, 51000, 71000];

% Pasar de geopotencial a geometrica
z_n_geom = (R_T .* z_n_geop) ./ (R_T - z_n_geop);

% Determinar la capa (n*)
n_star = find(z_star >= z_n_geom(1:end-1) & z_star < z_n_geom(2:end), 1, 'first');
 if isempty(n_star)
   if z_star < z_n_geom(1)
      n_star = 1;
   else
      n_star = numel(z_n_geom)-1;
   end
 end

% Definir los tramos del ascenso: z*_boundary(1) = 12000; z*_boundary(2) = z_n_geom(n_star+1); z*_boundary(3) = z_n_geom(n_star+2); ...
z_star_boundary = [z_star,z_n_geom(n_star+1:end)];
K = numel(z_star_boundary); % numero de tramos por encima de z_star

% Presión de referencia en cada frontera (P_prima = P^)
P_prima = zeros(1, K);

% Radio de referencia en cada frontera
R = zeros(1, K);

% Para el tramo 1 (k=1)
P_prima(1) = P(z_star);     % P^_0

Tz = T_globo(z_star);       % K
Pz = P_prima(1);             % Pa
R(1) = ( (3*m_he*R_eff)/(4*pi) * (Tz / Pz) )^(1/3);

% Resto de tramos
for k = 2:K
   P_prima(k) = P(z_star_boundary(k));
   R(k) = R(k-1) * ( P_prima(k-1) / P_prima(k) )^(2/5);
end

% Calcular R_globo(z) para cada altitud
R_g = zeros(size(z));

% Para z <= z*:
zona_baja = (z <= z_star);
 if any(zona_baja)
   T_z = T_globo(z(zona_baja));   
   P_z = P(z(zona_baja));         
   R_g(zona_baja) = ((3*m_he*R_eff)/(4*pi).*(T_z./P_z)).^(1/3);
 end

% Para z > z*
 zona_alta = (z > z_star);

    if any(zona_alta)
        z_alta = z(zona_alta);
        R_alta = zeros(size(z_alta));

        % Para cada punto de z_alta necesitamos averiguar en qué tramo está.
        % z*_boundary(k) <= z < z*_boundary(k+1)

        for i = 1:numel(z_alta)
            zi = z_alta(i);
            k = find( zi >= z_star_boundary, 1, 'last' );

            if isempty(k)
              k = 1;
            end

            % Constantes del tramo
            Pk = P_prima(k);   % P^_k
            Rk = R(k);             % R_k

            % Presión a esa altura
            P_actual = P(zi);

            R_alta(i) = Rk * ( Pk / P_actual )^(2/5);
        end

        R_g(zona_alta) = R_alta;
    end


end

