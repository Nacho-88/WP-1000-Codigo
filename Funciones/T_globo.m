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
 n_star = find(z_star >= z_n_geom(1:end-1) & z_star < z_n_geom(2:end), 1, 'first');
 if isempty(n_star)
    if z_star < z_n_geom(1)
      n_star = 1;
    else
      n_star = numel(z_n_geom)-1;
    end
 end

% Definir los tramos del ascenso: z*_frontera(1) = 12000; z*_frontera(2) = z_n_geom(n_star+1); z*_frontera(3) = z_n_geom(n_star+2); ...
z_star_frontera = [z_star,z_n_geom(n_star+1:end)];
K = numel(z_star_frontera); % numero de tramos por encima de z_star

% Presión de referencia en cada frontera (P_prima = P^)
P_prima = zeros(1, K);

% Temperatura de referencia en cada frontera (T_prima = T^)
T_prima = zeros(1, K);

% Para el tramo 1 (k=1)
P_prima(1) = P(z_star);     % P^_0
T_prima(1) = T_atm(z_star); % T^_0

% Resto de tramos
for k = 2:K
    P_prima(k) = P(z_star_frontera(k));
    T_prima(k) = T_prima(k-1) * ( P_prima(k) / P_prima(k-1) )^(2/5); % T^_(k) = T^_(k-1) * (P^_(k)/P^_(k-1))^(2/5)
end

% Calcular T_globo(z) para cada altitud
T = zeros(size(z));

% Para z <= z*: T_globo = T_atm(z)
zona_baja = (z <= z_star);

if any(zona_baja)
  T(zona_baja) = T_atm(z(zona_baja));
end

% Para z > z*
zona_alta = (z > z_star);

if any(zona_alta)
  z_alta = z(zona_alta); 
  T_alta = zeros(size(z_alta));

  % z*_frontera(k) <= z < z*_frontera(k+1)

  for i = 1:numel(z_alta)
     zi = z_alta(i);
     k = find( zi >= z_star_frontera, 1, 'last' );
     
     if isempty(k)
       k = 1;
     end

     % Constantes
     Pk = P_prima(k);     % P^_k
     Tk = T_prima(k);     % T^_k

     % Presión a esa altura
     P_actual = P(zi);

     % T_globo(zi) = T^_k * ( P(zi) / P^_k )^(2/5)
     T_alta(i) = Tk * ( P_actual / Pk )^(2/5);

   end

     % Guardar en vector
     T(zona_alta) = T_alta;

end



end

