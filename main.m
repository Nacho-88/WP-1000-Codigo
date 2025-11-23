%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOMBRES VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Altitud inicial: z0
% Altitud: z
% Gravedad a nivel del suelo: g_0
% Empuje: E
% Masa Helio: m_He
% Temperatura atmosferica: T_atm
% Presion atmosferica: P (Pa)
% Radio de la Tierra: R_T
% Rozamiento: F_roz
% Fuerza gravitatoria: F_g
% Constante de gravitación universal: G
% Masa de la Tierra: M_T
% Masa total: M
% Masa paracaídas: M_paraca
% Coeficiente de arrastre: C_D 
% Radio globo: R_globo
% Constante de gases ideales específica del aire: R
% Constante de gases ideales específica del Helio: R_prima
% Masa molecular promedio del aire: m_a
% Masa atómica del Helio: m_He_molar
% Distancia vertical entre globo-payload: l_caja_globo

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


M_caja = 3.2;      % Masa del payload (kg)

M_globo = 2;       % Masa del globo (kg)

M_paraca1 = 0.23;  % Masa paracaídas 1 (kg)

M_paraca2 = 0.08;  % Masa paracaídas 2 (kg)


save constantes G M_T R_T g_0 P_0 T_n R m_a m_He_molar R_prima z_n T_n C_D_caja C_D_globo C_D_paraca_1 C_D_paraca_2 l_caja_paraca1 l_paraca1_paraca2 l_paraca2_globo l_caja_globo L_x L_y L_z A_caja A1 A2 M_caja M_globo M_paraca1 M_paraca2 

%% Parámetros variables

R_exp = 12.4/2;                          % Radio de explosión (m)
z_exp_est = 37550;                       % Altitud de explosión estimada (m)
R_0 = 2.4/2;                             % Radio inicial (m)
T_amb = 18;                              % Temperatura ambiente (ºC)
P_amb = 89876;                           % Presión ambiental (Pa)
validacion = true;                       % booleano para elegir la forma de calcular la masa de helio


ruta = addpath('Funciones');

P_n = generar_P_n();      % Pa
save constantes.mat P_n -append


m_He = calcularMasaHelio(R_exp, R_0, z_exp_est, T_amb, P_amb, validacion);   % kg

[P_b, delta_P, Nozzle_Lift] = llenado(m_He, T_amb);  % P_b [bar]     Nozzle_Lift [N]


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


% Definición condiciones iniciales y parámetros para Runge-Kutta/ode45:
% dt = 1.05;    % s
t0 = 0;         % s
tf = 3*3600;    % s

% Aproximamos el radio del globo en altitud inicial a la altitud de la
% parte inferior del globo (pocos metro de diferencia y la presión a esa altitud varía poco)
z0 = L_z + l_caja_globo + R_globo(L_z+l_caja_globo, m_He);    % m

dz_dt0 = 0;     % m/s

w0 = [dz_dt0, z0+1026]';

f = @(t, w) ec_mov(t, w, m_He, R_exp);

% Uso Runge-Kutta/ode45 para la resolución del sistema de ecuaciones diferenciales
% [t, w] = Metodo_RK4(dt, t0, tf, f, w0);
[t, w] = ode45(f, [t0 tf], w0, odeset(RelTol=1e-5,AbsTol=1e-7));


% Desglose de la matriz devuelta por RK/ode45 en velocidades y altitudes (vectores fila)
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
save Resultados_simulaciones\vuelo_n11.mat z dz_dt t m_He z_exp t_exp  % Guardamos los parámetros de este vuelo en una carpeta dedicada a él.

% Vuelta a la ruta de búsqueda inicial (no la volvemos a usar en adelante).
path(ruta)
clear ruta tf_add t_aux w_aux index z_aux t0_aux



%% Gráficas y valores importantes
% =================================================================
% ANÁLISIS DE DATOS Y CÁLCULOS PRELIMINARES
% =================================================================
% Variables de entrada asumidas: 't' (tiempo), 'z' (altitud) y 'dz_dt' (velocidad).
z0 = z(1); % Altitud inicial del sistema.

% CÁLCULO DEL TIEMPO TOTAL DE VUELO
t_total = t(end); % Tiempo total del vuelo en segundos.
% Conversión del tiempo total a formato H:M:S para visualización.
H_total = floor(t_total / 3600);
M_total = floor(mod(t_total, 3600) / 60);
S_total = mod(t_total, 60);
t_str_total = sprintf('%02d h,%02d min, %05.2f s', H_total, M_total, S_total);

% =================================================================
% GRÁFICAS DE VUELO
% =================================================================

% --- Gráfica 1: Altitud vs. Tiempo (Perfil de Vuelo) ---
figure(1)
plot(t, z, 'b', 'MarkerSize', 8);
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

% Conversión del tiempo de explosión (t_max) a formato H:M:S para la etiqueta.
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
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.
hold off

% --- Gráfica 2: Velocidad vs. Tiempo (Perfil de Velocidades) ---
figure(2)
plot(t, dz_dt, 'b', 'MarkerSize', 8)
title('Perfil de Velocidades')
xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')
grid on
hold off

% =================================================================
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

% 1. Velocidad media de ascenso simple (Altura de explosión / Tiempo explosión).
v_avg_simple = (z_exp-z0) / (t_exp-t(1));
v_avg_suma_arit = sum(dz_dt(1:length(t_ascend)))/length(t_ascend);
v_avg_suma_pond = sum((dz_dt(2:length(t_ascend))+dz_dt(1:length(t_ascend)-1)).*(t_ascend(2:end)-t_ascend(1:end-1))./2)/(t_exp-t(1));

% 2. Velocidad media de ascenso por ajuste lineal (Mínimos Cuadrados).
[P_fit_1, S_fit_1, MU_fit_1] = polyfit(t_ascend, z_ascend, 1);
[z_fit_1, delta_fit_1] = polyval(P_fit_1, t_ascend, S_fit_1, MU_fit_1);
v_avg_fit_1 = P_fit_1(1)/MU_fit_1(2); % El coeficiente P(1) es la pendiente (velocidad). 
% Al calcular S y MU el polinomio devuelto hay que evaluarlo en (t-MU(1))/MU(2),
% por lo que la nueva pendiente es P(1)/MU(2).

% 3. Velocidad media de ascenso por ajuste parabólico.
[P_fit_2, S_fit_2, MU_fit_2] = polyfit(t_ascend, z_ascend, 2);
[z_fit_2, delta_fit_2] = polyval(P_fit_2, t_ascend, S_fit_2, MU_fit_2);
% Como z(t) = P_fit_2(1)*t^2 + P_fit_2(2)*t + P_fit_2(3) entonces:
% v(t) = 2*P_fit_2(1)*t + P_fit_2(2), así que el promedio queda
% <v> = P_fit_2(1)*(t_f + t_i) + P_fit_2(2)
% Hay que hacer la corrección debida a que al calcular S y MU el polinomio
% devuelto está evaluado en (t-MU(1))/MU(2) en vez de t.
v_avg_fit_2 = P_fit_2(1)/(MU_fit_2(2))^2*(t_ascend(end)+t_ascend(1)-2*MU_fit_2(1)) + P_fit_2(2)/MU_fit_2(2);

v_inicial_min_cuad = P_fit_2(2)/MU_fit_2(2) - 2*P_fit_2(1)*MU_fit_2(1)/(MU_fit_2(2))^2;
a_min_cuad = P_fit_2(1)/(2*MU_fit_2(2)^2);

% Tiempo de transición
t_trans = t_max - t_exp;

% Conversión a formato H:M:S
H_trans = floor(t_trans / 3600);
M_trans = floor(mod(t_trans, 3600) / 60);
S_trans = mod(t_trans, 60);
t_str_trans = sprintf('%02d h,%02d min, %05.2f s', H_trans, M_trans, S_trans);

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
fprintf('Duración del periodo de transición: %s\n', t_str_trans);
fprintf('-------------------------------------------\n');
fprintf('Velocidad media de ascenso (simple): %.4f m/s\n', v_avg_simple);
fprintf('Velocidad media de ascenso (media aritmética): %.4f m/s\n', v_avg_suma_arit);
fprintf('Velocidad media de ascenso (media ponderada): %.4f m/s\n', v_avg_suma_pond);
fprintf('Velocidad media de ascenso (ajuste lineal): %.4f m/s\n', v_avg_fit_1);
fprintf('Valor de R\xB2 para el ajuste lineal: %.6f \n', S_fit_1.rsquared);
fprintf('Velocidad media de ascenso (ajuste parabólico): %.4f m/s\n', v_avg_fit_2);
fprintf('Valor de R\xB2 para el ajuste parabólico: %.6f \n', S_fit_2.rsquared);

% Gráficas de los ajustes de mínimos cuadrados para el ascenso
figure(3)
hold on
plot(t_ascend, z_ascend, 'b', 'MarkerSize', 8);
plot(t_ascend, z_fit_1, 'g--', 'LineWidth', 2);
plot(t_ascend, z_fit_2, 'm--', 'LineWidth', 2);
hold off

str_fit_1 = sprintf('Ajuste lineal (R\xB2 = %.4f)', S_fit_1.rsquared);
str_fit_2 = sprintf('Ajuste parabólico (R\xB2 = %.4f)', S_fit_2.rsquared);

title('Altitud durante el Ascenso con Ajustes Polinomiales')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Ascenso', str_fit_1, str_fit_2, 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


figure(4)
hold on
plot(t_ascend, z_ascend, 'b', 'MarkerSize', 8);
plot(t_ascend, z_fit_1, 'g', 'LineWidth', 2);
plot(t_ascend, z_fit_1 + delta_fit_1, 'm--', 'LineWidth', 2);
plot(t_ascend, z_fit_1 - delta_fit_1, 'm--', 'LineWidth', 2);
plot(t_ascend, z_fit_1 + 2*delta_fit_1, 'r--', 'LineWidth', 2);
plot(t_ascend, z_fit_1 - 2*delta_fit_1, 'r--', 'LineWidth', 2);
hold off

title('Ajuste polinomial (Grado 1) durante el ascenso')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Ascenso', str_fit_1, 'Margen de error de 68%', '', 'Margen de error de 95%', '', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


figure(5)
hold on
plot(t_ascend, z_ascend, 'b', 'MarkerSize', 8);
plot(t_ascend, z_fit_2, 'g', 'LineWidth', 2);
plot(t_ascend, z_fit_2 + delta_fit_2, 'm--', 'LineWidth', 2);
plot(t_ascend, z_fit_2 - delta_fit_2, 'm--', 'LineWidth', 2);
plot(t_ascend, z_fit_2 + 2*delta_fit_2, 'r--', 'LineWidth', 2);
plot(t_ascend, z_fit_2 - 2*delta_fit_2, 'r--', 'LineWidth', 2);
hold off

title('Ajuste polinomial (Grado 2) durante el ascenso')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Ascenso', str_fit_2, 'Margen de error de 68%', '', 'Margen de error de 95%', '', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


% =================================================================
% TAREA 3: MODELADO DEL DESCENSO CON AJUSTES FUNCIONALES
% =================================================================

% Datos de la fase de descenso (desde t_max hasta el final).
t_descent = t(idx_max:end);
z_descent = z(idx_max:end);
% Normalizar el tiempo para que t'=0 en el punto máximo.
t_prime_descent = t_descent - t_descent(1);

% 1. Ajuste Polinomial (Grado 2).
[P_poly2, S_2, MU_poly_2] = polyfit(t_prime_descent, z_descent, 2);
% Evaluar el polinomio usando los parámetros de normalización (MU) para mayor estabilidad.
[z_fit_poly2, delta_2] = polyval(P_poly2, t_prime_descent, S_2, MU_poly_2);

% 2. Ajuste Polinomial (Grado 3).
[P_poly3, S_3, MU_poly_3] = polyfit(t_prime_descent, z_descent, 3);
% Evaluar el polinomio usando los parámetros de normalización (MU) para mayor estabilidad.
[z_fit_poly3, delta_3] = polyval(P_poly3, t_prime_descent, S_3, MU_poly_3);

% 3. Ajuste Exponencial (Aproximación lineal sobre el logaritmo).
z_final = z(end);
z_diff = z_descent - z_final;
% Evitar log(0): reemplazar valores no positivos con un número muy pequeño (eps).
z_diff(z_diff <= 0) = eps;

% Ajuste lineal (Grado 1) sobre log(z - z_final) vs. t'.
[P_exponent, S_exponent, MU_exponent] = polyfit(t_prime_descent, log(z_diff), 1);
[z_fit_exponent, delta_exponent] = polyval(P_exponent, t_prime_descent, S_exponent, MU_exponent);


% --- Gráficas 6, 7, 8 y 9: Altitud en el Descenso con Ajustes ---
figure(6)
hold on
plot(t_descent, z_descent, 'b', 'MarkerSize', 8);
plot(t_descent, z_fit_poly2, 'g--', 'LineWidth', 2);
plot(t_descent, z_fit_poly3, 'm--', 'LineWidth', 2);
plot(t_descent, (z_final + exp(z_fit_exponent)), 'r--', 'LineWidth', 2);
hold off

str_fit_poly2 = sprintf('Ajuste parabólico (R\xB2 = %.4f)', S_2.rsquared);
str_fit_poly3 = sprintf('Ajuste cúbico (R\xB2 = %.4f)', S_3.rsquared);
str_fit_exponent = sprintf('Ajuste exponencial (R\xB2 = %.4f)', S_exponent.rsquared);

title('Altitud durante el Descenso con Ajustes (Polinomiales y Exponencial)')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Descenso', str_fit_poly2, str_fit_poly3, str_fit_exponent, 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


figure(7)
hold on
plot(t_descent, z_descent, 'b', 'MarkerSize', 8);
plot(t_descent, z_fit_poly2, 'g', 'LineWidth', 2);
plot(t_descent, z_fit_poly2 + delta_2, 'm--', 'LineWidth', 2);
plot(t_descent, z_fit_poly2 - delta_2, 'm--', 'LineWidth', 2);
plot(t_descent, z_fit_poly2 + 2*delta_2, 'r--', 'LineWidth', 2);
plot(t_descent, z_fit_poly2 - 2*delta_2, 'r--', 'LineWidth', 2);
hold off

title('Ajuste polinomial (grado 2) durante el descenso')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Descenso', str_fit_poly2, 'Margen de error de 68%', '', 'Margen de error de 95%', '', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


figure(8)
hold on
plot(t_descent, z_descent, 'b', 'MarkerSize', 8);
plot(t_descent, z_fit_poly3, 'g', 'LineWidth', 2);
plot(t_descent, z_fit_poly3 + delta_3, 'm--', 'LineWidth', 2);
plot(t_descent, z_fit_poly3 - delta_3, 'm--', 'LineWidth', 2);
plot(t_descent, z_fit_poly3 + 2*delta_3, 'r--', 'LineWidth', 2);
plot(t_descent, z_fit_poly3 - 2*delta_3, 'r--', 'LineWidth', 2);
hold off

title('Ajuste polinomial (grado 3) durante el descenso')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Descenso', str_fit_poly3, 'Margen de error de 68%', '', 'Margen de error de 95%', '', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


figure(9)
hold on
plot(t_descent, z_descent, 'b', 'MarkerSize', 8);
plot(t_descent, (z_final + exp(z_fit_exponent)), 'g', 'LineWidth', 2);
plot(t_descent, (z_final + exp(z_fit_exponent + delta_exponent)), 'm--', 'LineWidth', 2);
plot(t_descent, (z_final + exp(z_fit_exponent - delta_exponent)), 'm--', 'LineWidth', 2);
plot(t_descent, (z_final + exp(z_fit_exponent + 2*delta_exponent)), 'r--', 'LineWidth', 2);
plot(t_descent, (z_final + exp(z_fit_exponent - 2*delta_exponent)), 'r--', 'LineWidth', 2);
hold off

title('Ajuste exponencial durante el descenso')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos de Descenso', str_fit_exponent, 'Margen de error de 68%', '', 'Margen de error de 95%', '', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.



% Gráfica 10: comparación final entre curva simulada y curvas de mínimos
% cuadrados consideradas
z_min_cuad = [z_fit_2; z_fit_poly3];
z_asc_des = [z_ascend; z_descent];
t_asc_des = [t_ascend; t_descent];

figure(10)
hold on
plot(t_asc_des, z_asc_des, 'b', 'MarkerSize', 8);
plot(t_asc_des, z_min_cuad, 'r--', 'LineWidth', 2);
hold off

title('Comparación datos con ajustes')
xlabel('Tiempo (s)')
ylabel('Altitud (m)')
legend('Datos Ascenso-Descenso', 'Ajustes por Mínimos Cuadrados', 'Location', 'best')
grid on
ax = gca;
ax.YAxis.Exponent = 0; % Evitar notación científica en el eje Y.


%% Limpieza variables intermedias

clear H H_fall H_total H_trans
clear M M_fall M_total M_trans
clear S S_fall S_total S_trans
clear t_ascend t_asc_des t_descent t_exp_str t_fall_duration t_max t_prime_descent t_str t_str_fall t_str_total t_str_trans t_trans t_total
clear delta_fit_1 MU_fit_1 P_fit_1 S_fit_1 z_fit_1 str_fit_1
clear delta_fit_2 MU_fit_2 P_fit_2 S_fit_2 z_fit_2 str_fit_2
clear delta_2 MU_poly_2 P_poly2 S_2 z_fit_poly2 str_fit_poly2
clear delta_3 MU_poly_3 P_poly3 S_3 z_fit_poly3 str_fit_poly3
clear delta_exponent MU_exponent P_exponent S_exponent z_fit_exponent str_fit_exponent
clear z_ascend z_asc_des z_descent z_diff z_min_cuad z_final
clear ax idx_max