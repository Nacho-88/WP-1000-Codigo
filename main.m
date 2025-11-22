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
%Tglobo 1 es igual a la atm 
%radio globo= r_globo

%% Constantes
G   = 6.67e-11; % N * m^2 / (kg^2)
M_T = 5.97e24;  % kg
R_T = 6371000;  % Radio de la Tierra medio (m)
g_0 = 9.80665;  % m / s^2
P_0 = 101325;   % Pa


R = 287.053;                       % J/(kg·K)

m_a = 0.0289644;                   % kg/mol

m_He_molar = 0.0040026;            % kg/mol

R_prima = (m_a / m_He_molar) * R;  % J/(kg*K)


% Separación en capas del modelo
% La última capa se extiende artificialmente para que no haya errores en los pasos intermedios del ode45
z_n = [0, 11000, 20000, 32000, 47000, 51000, 71000, 170000];    % Altitudes (geopotenciales) (m)

T_n = [15, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5, -256,5];     % Temperaturas (ºC)
T_n = T_n + 273.15*ones(size(T_n));                     % Temperaturas (K)


% Coeficientes de arrastre (adimensionales)
C_D_caja = 1.57439;   % Coeficiente de arrastre caja

C_D_globo = 0.32;      % Coeficiente de arrastre globo

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


M_caja = 2.21553;  % Masa del payload (kg)

M_globo = 2;       % Masa del globo (kg)

M_paraca1 = 0.23;  % Masa paracaídas 1 (kg)

M_paraca2 = 0.08;  % Masa paracaídas 2 (kg)


save constantes G M_T R_T g_0 P_0 T_n R m_a m_He_molar R_prima z_n T_n C_D_caja C_D_globo C_D_paraca_1 C_D_paraca_2 l_caja_paraca1 l_paraca1_paraca2 l_paraca2_globo l_caja_globo L_x L_y L_z A_caja A1 A2 M_caja M_globo M_paraca1 M_paraca2 

%% Parámetros variables

R_exp = 12.4/2;                          % m
z_exp_est = 37550;                       % m
R_0 = 2.55/2;                            % m
T_amb = 15;                              % ºC
P_amb = 89876;                           % Pa
validacion = true;                       % booleano para elegir la forma de calcular la masa de helio


ruta = addpath('Funciones');

P_n = generar_P_n();      % Pa
save constantes.mat P_n -append


m_He = calcularMasaHelio(R_exp, R_0, z_exp_est, T_amb, P_amb, validacion);   % kg

[P_b, delta_P, Nozzle_Lift] = llenado(m_He, R_0, T_amb);  % P_b [bar]     Nozzle_Lift [N]


path(ruta);


%% Resolución ecuación del movimiento

% Cambio de la ruta de búsqueda para poder buscar en la carpeta de funciones.
ruta = addpath('Funciones');

% En caso de haber ejecutado previamente el código, tenemos que eliminar
% las variables t_exp y z_exp del archivo parametros.mat para poder
% detectar correctamente el instante de explosión en la próxima simulación.
S = load('parametros.mat');
if isfield(S,'t_exp') | isfield(S, 'z_exp')
    S = rmfield(S, "z_exp");
    S = rmfield(S, "t_exp");
    save("parametros.mat", "-struct", "S")
end
clear S


% Definición condiciones iniciales y parámetros para Runge-Kutta:
% dt = 1.05;    % s
t0 = 0;         % s
tf = 3*3600;    % s

% Aproximamos el radio del globo en altitud inicial a la altitud de la
% parte inferior del globo (pocos metro de diferencia y la presión a esa altitud varía poco)
z0 = L_z + l_caja_globo + R_globo(L_z+l_caja_globo, m_He);    % m

dz_dt0 = 0;     % m/s

w0 = [dz_dt0, z0+1026]';

f = @(t, w) ec_mov(t, w, m_He, R_exp);

% Uso Runge-Kutta para la resolución del sistema de ecuaciones diferenciales
% [t, w] = Metodo_RK4(dt, t0, tf, f, w0);
[t, w] = ode45(f, [t0 tf], w0, odeset(RelTol=1e-5,AbsTol=1e-7));


% Desglose de la matriz devuelta por RK en velocidades y altitudes (vectores fila)
dz_dt = w(:,1); % m/s
z = w(:,2);     % m


% Si no hemos alcanzado altitud 0 (o próxima a 0) continuamos un tiempo
% extra hasta alcanzarlo:
while z(end)>(L_z/2)

    w0_aux = [dz_dt(end), z(end)]';
    t0_aux = tf;
    tf = tf + 900; % Tiempo añadido (en segundos) para alcanzar altitud 0
    [t_aux, w_aux] = ode45(f, [t0_aux tf], w0_aux, odeset(RelTol=1e-7,AbsTol=1e-9));

    % Añadimos los nuevos valores de tiempos y el vector de estados.
    t = [t; t_aux(2:end)];

    % Desglose de la matriz
    dz_dt = [dz_dt; w_aux(:,1)];
    z = [z; w_aux(:,2)];

    w = [dz_dt, z];

end


% Eliminamos los elementos con una altura negativa
z = z(1:length(z(z>=0))+1);
% Interpolamos linealmente para saber el instante de tiempo en el que se
% llega a z=0.
t(length(z)) = t(length(z)-1) + (t(length(z))-t(length(z)-1))*(0-z(end-1))/(z(end)-z(end-1));
z(end) = 0;
dz_dt(length(z)+1:end) = [];
t(length(z)+1:end) = [];

% Desplazamiento de la altitud para que corresponda al payload en vez del globo:
load parametros.mat z_exp t_exp
t_aux = t(1);
index = 1;
while (t_aux<t_exp)
    z(index) = z(index) - R_globo(z(index), m_He) - l_caja_globo - L_z/2;
    index = index + 1;
    t_aux = t(index);
end

save parametros z dz_dt t -append   % Guardamos los parámetros calculados.

% Hay que cambiar el archivo en función del vuelo que se esté simulando
save Funciones\vuelo_n8.mat z dz_dt t m_He z_exp t_exp  % Guardamos los parámetros de este vuelo en una carpeta dedicada a él.

% Vuelta a la ruta de búsqueda inicial (no la volvemos a usar en adelante).
path(ruta)
clear ruta tf_add t_aux w_aux index z_aux t0_aux



%% =================================================================
% ANÁLISIS DE DATOS Y CÁLCULOS PRELIMINARES
% =================================================================
% Variables de entrada asumidas: 't' (tiempo), 'z' (altitud) y 'dz_dt' (velocidad).
z0 = z(1); % Altitud inicial del sistema.

%% CÁLCULO DEL TIEMPO TOTAL DE VUELO
t_total = t(end); % Tiempo total del vuelo en segundos.
% Conversión del tiempo total a formato H:M:S para visualización.
H_total = floor(t_total / 3600);
M_total = floor(mod(t_total, 3600) / 60);
S_total = mod(t_total, 60);
t_str_total = sprintf('%02d h,%02d min, %05.2f s', H_total, M_total, S_total);

%% =================================================================
% GRÁFICAS DE VUELO
% =================================================================

% --- Gráfica 1: Altitud vs. Tiempo (Perfil de Vuelo) ---
figure(1)
plot(t, z);
hold on
grid on

% Cálculo del punto máximo (Altitud y Tiempo).
[z_max, idx_max] = max(z);
t_max = t(idx_max);

% Conversión del tiempo máximo (t_max) a formato H:M:S para la etiqueta.
H = floor(t_max / 3600);
M = floor(mod(t_max, 3600) / 60);
S = mod(t_max, 60);
t_str = sprintf('%02d h,%02d min, %05.2f s', H, M, S);

% Conversión del tiempo de explosiòn (t_max) a formato H:M:S para la etiqueta.
H = floor(t_exp / 3600);
M = floor(mod(t_exp, 3600) / 60);
S = mod(t_exp, 60);
t_exp_str = sprintf('%02d h,%02d min, %05.2f s', H, M, S);

% Marcar el punto de altitud máxima y añadir una etiqueta informativa.
plot(t_max, z_max, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
text(t_max, z_max, ...
    ['  Máximo (', t_str, ' y ', num2str(z_max, '%.0f'), ' m)'], ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'FontSize', 10);

title('Perfil de Ascenso y Descenso')
xlabel('tiempo (s)')
ylabel('altura (m)')
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.
hold off

% --- Gráfica 2: Velocidad vs. Tiempo (Perfil de Velocidades) ---
figure(2)
plot(t, dz_dt)
title('Perfil de Velocidades')
xlabel('tiempo (s)')
ylabel('velocidad (m/s)')
grid on
hold off

%% =================================================================
% CÁLCULOS ADICIONALES (Duración de Caída y Velocidades Medias)
% =================================================================

% --- Tarea 1: Duración de la caída del payload ---
% CÁLCULO DE LA DURACIÓN DE LA CAÍDA: Tiempo de recuperación - Tiempo en el máximo.
t_fall_duration = t_total - t_max;

% Conversión de la duración de caída a formato H:M:S.
H_fall = floor(t_fall_duration / 3600);
M_fall = floor(mod(t_fall_duration, 3600) / 60);
S_fall = mod(t_fall_duration, 60);
t_str_fall = sprintf('%02d h,%02d min, %05.2f s', H_fall, M_fall, S_fall);

% --- Tarea 2: Velocidades medias de ascenso ---
t_ascend = t(t<=t_exp); % Datos de tiempo durante el ascenso.
z_ascend = z(1:length(t_ascend)); % Datos de altitud durante el ascenso.

% 1. Velocidad media de ascenso simple (Altura Máxima / Tiempo Máximo).
v_avg_simple = (z_exp-z0) / t_exp;

% 2. Velocidad media de ascenso por ajuste lineal (Mínimos Cuadrados).
P_fit = polyfit(t_ascend, z_ascend, 1);
v_avg_fit = P_fit(1); % El coeficiente P(1) es la pendiente (velocidad).

% IMPRESIÓN DE RESULTADOS
fprintf('\n--- Resultados de Datos y Velocidad ---\n');
fprintf('Altitud inicial: %.2f m\n', z0);
fprintf('Tiempo de explosión del globo: %s\n', t_exp_str);
fprintf('Altitud de explosión: %.2f m\n', z_exp);
fprintf('Tiempo de máxima altitud del globo: %s\n', t_str);
fprintf('Altitud máxima alcanzada: %.2f m\n', z_max);
fprintf('-------------------------------------------\n');
fprintf('Tiempo total del vuelo (final de los datos): %s\n', t_str_total);
fprintf('Duración de la caída del payload: %s\n', t_str_fall);
fprintf('-------------------------------------------\n');
fprintf('Velocidad media de ascenso (simple): %.4f m/s\n', v_avg_simple);
fprintf('Velocidad media de ascenso (ajuste lineal): %.4f m/s\n', v_avg_fit);


%% =================================================================
% TAREA 3: MODELADO DEL DESCENSO CON AJUSTES FUNCIONALES
% =================================================================

% Datos de la fase de descenso (desde t_max hasta el final).
t_descent = t(idx_max:end);
z_descent = z(idx_max:end);
% Normalizar el tiempo para que t'=0 en el punto máximo.
t_prime_descent = t_descent - t_descent(1);

% 1. Ajuste Polinomial (Grado 3).
[P_poly3, ~, MU_poly] = polyfit(t_prime_descent, z_descent, 3);
% Evaluar el polinomio usando los parámetros de normalización (MU) para mayor estabilidad.
z_fit_poly3 = polyval(P_poly3, t_prime_descent, [], MU_poly);

% 2. Ajuste Exponencial (Aproximación lineal sobre el logaritmo).
z_final = z(end);
z_diff = z_descent - z_final;
% Evitar log(0): reemplazar valores no positivos con un número muy pequeño (eps).
z_diff(z_diff <= 0) = eps;

% Ajuste lineal (Grado 1) sobre log(z - z_final) vs. t'.
[P_exp, ~, MU_exp] = polyfit(t_prime_descent, log(z_diff), 1);

% Calcular los coeficientes c1 y c2 para el modelo exponencial z = z_final + exp(c1 + c2 * t').
c2_exp = P_exp(1) / MU_exp(2);      % c2 (pendiente)
c1_exp = P_exp(2) - c2_exp * MU_exp(1); % c1 (intercepto)

% Evaluar el modelo exponencial.
z_fit_exp = z_final + exp(c1_exp + c2_exp * t_prime_descent);

% --- Gráfica 3: Altitud en el Descenso con Ajustes ---
figure(3)
plot(t_descent, z_descent, 'b.', 'MarkerSize', 8);
hold on
plot(t_descent, z_fit_poly3, 'r-', 'LineWidth', 2);
plot(t_descent, z_fit_exp, 'g--', 'LineWidth', 2);
hold off

title('Altitud durante el Descenso con Ajustes (Polinomial y Exponencial)')
xlabel('tiempo (s)')
ylabel('altura (m)')
legend('Datos de Descenso', 'Ajuste Polinomial (Grado 3)', 'Ajuste Exponencial', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


