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
g_0 = 9.80665;  % m / s^2
P_0 = 101325;   % Pa


z_star = 12000;                    % Límite para cambio de formula (m)

R = 287.053;                       % J/(kg·K)

m_a = 0.0289644;                   % kg/mol

m_He_molar = 0.0040026;            % kg/mol

R_prima = (m_a / m_He_molar) * R;  % J/(kg*K)


% Separación en capas del modelo
z_n = [0, 11000, 20000, 32000, 47000, 51000, 71000];    % Altitudes (geopotenciales) (m)

T_n = [15, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5];     % Temperaturas (ºC)
T_n = T_n + 273.15*ones(size(T_n));                     % Temperaturas (K)


% Coeficientes de arrastre (adimensionales)
C_D_caja = 1.57439;   % Coeficiente de arrastre caja

C_D_globo = 0.3;      % Coeficiente de arrastre globo

C_D_paraca_1 = 0.97;  % Coeficiente de paracaidas 1

C_D_paraca_2 = 1.17;  % Coeficiente de paracaidas 2


% Distancias verticales (m)
l_caja_paraca1  = 12.0;    % Caja – Paracaídas 1

l_paraca1_paraca2 = 15.0;  % Paracaídas 1 – Paracaídas 2

l_paraca2_globo = 10.0;    % Paracaídas 2 – Globo

l_caja_globo = l_caja_paraca1 + l_paraca1_paraca2 + l_paraca2_globo;     % Longitud (m) de la cuerda entre globo y payload


% Dimensiones caja
L_x = 0.35; % (m)
L_y = 0.28; % (m)
L_z = 0.28; % (m)

A_caja = L_x * L_y;    % Área payload (m^2)

A1 = 3.25;         % Área paracaidas 1 (m^2)

A2 = 1.13;         % Área paracaidas 1 (m^2)


M_caja = 2;        % Masa del payload (kg)

M_globo = 2;       % Masa del globo (kg)

M_paraca1 = 0.23;  % Masa paracaídas 1 (kg)

M_paraca2 = 0.08;  % Masa paracaídas 2 (kg)


save constantes G M_T R_T g_0 P_0 T_n z_star R m_a m_He_molar R_prima z_n T_n C_D_caja C_D_globo C_D_paraca_1 C_D_paraca_2 l_caja_paraca1 l_paraca1_paraca2 l_paraca2_globo l_caja_globo L_x L_y L_z A_caja A1 A2 M_caja M_globo M_paraca1 M_paraca2 

%% Parámetros variables

R_exp = 12.4/2;                          % m
z_exp_est = 37550;                  % m


ruta = addpath('Funciones');

P_n = generar_P_n();      % Pa
save constantes.mat P_n -append

% P_k = generar_P_k();	  % Pa
% save constantes.mat P_k -append
% 
% T_k = generar_T_k();      % K
% save constantes.mat T_k -append

m_He = calcularMasaHelio(R_exp, z_exp_est);   % kg

% R_k = generar_R_k(m_He);  % m
% save constantes.mat R_k m_He -append

path(ruta);


%% Resolución ecuación del movimiento

% Cambio de la ruta de búsqueda para poder buscar en la carpeta de funciones.
ruta = addpath('Funciones');

% Dado que se empieza con velocidad 0, el rozamiento no se tiene en cuenta
% en el primer cálculo, esto da una aceleración muy elevada que debería
% durar poco tiempo, por lo que vamos a considerar un paso pequeño para los
% primeros 15 puntos y después aumentaremos el paso ya que las
% aceleraciones serán pequeñas.

% Definición condiciones iniciales y parámetros para Runge-Kutta:
dt = 1.05;         % s
t0 = 0;         % s
tf = 4*3600;    % s

% Aproximamos el radio del globo en altitud inicial a la altitud de la
% parte inferior del globo (pocos metro de diferencia y la presión a esa altitud varía poco)
z0 = L_z + l_caja_globo + R_globo(L_z+l_caja_globo, m_He);    % m

dz_dt0 = 0;     % m/s

w0 = [dz_dt0, z0]';

f = @(t, w) ec_mov(t, w, m_He, R_exp);

% Uso Runge-Kutta para la resolución del sistema de ecuaciones diferenciales
% [t, w] = Metodo_RK4(dt, t0, tf, f, w0);
[t, w] = ode45(f, [t0 tf], w0, odeset(RelTol=1e-5,AbsTol=1e-7));


% Desglose de la matriz devuelta por RK en velocidades y altitudes (vectores fila)
dz_dt = w(:,1); % m/s
z = w(:,2);     % m


% Si no hemos alcanzado altitud 0 (o próxima a 0) continuamos un tiempo
% extra hasta alcanzarlo:
if z(end)>(L_z/2)
    w0_aux = [dz_dt(end), z(end)]';
    tf_add = tf + 3600; % Tiempo añadido (en segundos) para alcanzar altitud 0
    [t_aux, w_aux] = ode45(f, [tf tf_add], w0_aux, odeset(RelTol=1e-5,AbsTol=1e-7));
    
    % Añadimos los nuevos valores de tiempos y el vector de estados.
    t = [t, t_aux(2:end)];
    w = [w(:,1), w_aux(:,1); w(:,2), w_aux(:,2)];

    % Desglose de la matriz
    dz_dt = w(:,1);
    z = w(:,2);
end


% Desplazamiento de la altitud para que corresponda al payload en vez del globo:
load parametros.mat z_exp t_exp
z_aux = z(1);
index = 1;
while (z_aux<z_exp)
    z(index) = z(index) - R_globo(z(index), m_He) - l_caja_globo - L_z/2;
    index = index + 1;
    z_aux = z(index);
end

save parametros z dz_dt t -append   % Guardamos los parámetros calculados.


% Vuelta a la ruta de búsqueda inicial (no la volvemos a usar en adelante).
path(ruta)
clear ruta tf_add t_aux w_aux index z_aux



%% Gráficas y datos importantes a comparar
figure(1)
plot(t, z);
title('Perfil de ascenso y descenso')
xlabel('tiempo (s)')
ylabel('altura (m)')

figure(2)
plot(t,dz_dt)
title('Perfil de velocidades')
xlabel('tiempo (s)')
ylabel('velocidad (m/s)')

