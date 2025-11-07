%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOMBRES VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%altura inicial z0
%altura z
%gravedad g_0
%Empuje E
%masa helio m_He
%temperatura atmosferica T_atm
%temperatura globo T_globo
%presion atmosferica P (Pa)
%Radio tierra R_T
%rozamiento F_r
%G
%M_T
%masa total M
%masa paracaidas m_paracaidas
%coeficiente de arrastre c_d=0.75 d=1.64m (paracaidas)
%Tglobo 1 es igual a la atm hasta 12 km y luego Tglobo 2 varía a partir de
%esa altitud
%radio globo= r_globo

%% Constantes
G   = 6.67e-11; % N * m^2 / (kg^2)
M_T = 5.97e24;  % kg
R_T = 6371000;  % Radio de la Tierra medio (m)
g_0 = 9.81;     % m / s^2
P_0 = 101325;   % Pa

z_star = 12000; % Límite para cambio de formula (m)

R = 287.053;  % J/(kg·K)

m_a = 0.0289644;        % kg/mol

m_He_molar = 0.0040026;   % kg/mol


% T_n = [15, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5]; % ºC
% T_n = T_n + 273.15.*ones(size(T_n));        % K
% z_n = 1e3 .* [0, 11, 20, 32, 47, 51, 71];   % m

z_n_geop = [0, 11000, 20000, 32000, 47000, 51000, 71000];   % Vectores Altitud (geopotencial)

T_n = [288, 216.5, 216.5, 228.5, 270.5, 270.5, 214.5];      % Temperatura (K)

z_n_geom = (R_T .* z_n_geop) ./ (R_T - z_n_geop);           % Pasar de geopotencial a geometrica


R_exp = 12; % m

m_He = 2;   % kg

% r_globo=1;
% c_d=0.75;
dpar=1.64;
m_atomicaHe=4;
% Cd=0.3;
v=5.33;

C_D_caja = 1.57439;   % Coeficiente de arrastre caja

C_D_globo = 0.3;      % Coeficiente de arrastre globo

C_D_paraca_1 = 0.97;  % Coeficiente de paracaidas 1

C_D_paraca_2 = 1.17;  % Coeficiente de paracaidas 2


l_caja_globo = 37;         % Longitud (m) de la cuerda entre globo y payload

l_caja_paraca1  = 12.0;    % Caja – Paracaídas 1

l_paraca1_paraca2 = 15.0;  % Paracaídas 1 – Paracaídas 2

l_paraca2_globo = 10.0;    % Paracaídas 2 – Globo


A_caja = 0.098;    % Area payload [m^2]

A1 = 3.25;         % Area paracaidas 1

A2 = 1.13;         % Area paracaidas 1


M_caja = 2;        % Masa del payload (kg)

M_globo = 2;       % Masa del globo (kg)

M_paraca1 = 0.23;  % Masa paracaidas 1

M_paraca2 = 0.08;   % Masa paracaidas 1

% %Fórmula temperatura atmosférica
% Tp = @(z) Tn + (z - zn) * (Tnp1 - Tn) / (znp1 - zn);
% 
% altitud_m = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40] * 1000;
% temperatura_K = [15, 2, -10.98, -23.96, -36.93, -49.9, -56.5, -56.5, -56.5, -56.5, ...
%                  -56.5, -54.5, -52.5, -50.5, -48.62, -46.64, -44.66, -39.41, -33.87, -28.33, -22.8] + 273.15;

save constantes G R_T M_T g_0 P_0 T_n z_n



%% Resolución ecuación del movimiento

% Cambio de la ruta de búsqueda para poder buscar en la carpeta de funciones.
ruta = addpath('Funciones');


% Definición condiciones iniciales y parámetros para Runge-Kutta:
dt = 60;        % s
t0 = 0;         % s
tf = 3*3600;    % s
z0 = 0;         % m
dz_dt0 = 0;     % m/s

w0 = [dz_dt0, z0]';

f = @(t, w) ec_mov(t, w, m_He, R_exp);

% Uso Runge-Kutta para la resolución del sistema de ecuaciones diferenciales
[t, w] = Metodo_RK4(dt, t0, tf, f, w0);

% Desglose de la matriz devuelta por RK en velocidades y altitudes (vectores fila)
dz_dt = w(1,:); % m/s
z = w(2,:);     % m


% Si no hemos alcanzado altitud 0 (o próxima a 0) continuamos un tiempo
% extra hasta alcanzarlo:
if z(end)>(L_z/2)
    tf_add = 3600; % Tiempo añadido (en segundos) para alcanzar altitud 0
    [t_aux, w_aux] = Metodo_RK4(dt, tf, tf+tf_add, f, w(:,end));
    
    % Añadimos los nuevos valores de tiempos y el vector de estados.
    t = [t, t_aux(2:end)];
    w = [w(1,:), w_aux(1,:); w(2,:), w_aux(2,:)];

    % Desglose de la matriz
    dz_dt = w(1,:);
    z = w(2,:);
end


% Desplazamiento de la altitud para que corresponda al payload en vez del globo:
load parametros.mat z_exp t_exp
z_aux = z(1);
index = 1;
while (z_aux<z_exp)
    z(index) = z(index) - R_globo(z(index)) - l_c;
    index = index + 1;
    z_aux = z(index);
end

save parametros z dz_dt t   % Guardamos los parámetros calculados.


% Vuelta a la ruta de búsqueda inicial (no la volvemos a usar en adelante).
path(ruta)
clear ruta tf_add t_aux w_aux index z_aux



%% Gráficas y datos importantes a comparar


