%% Descripción
% Queremos hacer una estimación de la velocidad de ascenso a partir de 
% promedios. Para ello, aproximamos que la aceleración de ascenso se 
% mantiene constante y, por tanto, es igual a su promedio. A partir de este
% valor de la aceleración se deduce la velocidad promedio.


%% Parámetros iniciales
z_0 = 1026;     % Altitud inicial (m)
dz = 10;        % Paso de la integración numérica (m)

R_0 = 2.25/2;   % Radio inicial del globo (m)
z_exp = 37550;  % Altitud de explosión del globo (m)

validacion = true;

T_amb = 18;     % ºC
P_amb = 89876;  % Pa

ruta = addpath('Funciones');
m_He = calcularMasaHelio(NaN, R_0, NaN, T_amb, P_amb, validacion);
path(ruta)


%% Cálculo directo del promedio de la aceleración hasta explosión
% Promediando (respecto a las altitudes) la ecuación de movimiento del
% globo y teniendo en cuenta la aproximación de la aceleración por su
% promedio (a=cte => (dz_dt)^2=2*a*(z-z_0)), podemos simplificar la ecuación
% y despejar directamente la aceleración media.

z_f = z_exp;    % El promedio se hace desde z=z_0 hasta z=z_f

ruta = addpath('Funciones');
load constantes.mat G M_T R_T M_caja l_caja_globo L_z C_D_globo C_D_caja L_x L_y M_paraca1 M_paraca2 M_globo
F_g_caja_z = @(z) G*M_T*M_caja/((R_T+z-R_globo(z,m_He)-l_caja_globo-L_z)*(R_T+z-R_globo(z,m_He)));
F_roz_norm_z = @(z) -(z-z_0)*(pi*C_D_globo*(R_globo(z,m_He))^2*rho_atm(z+3*R_globo(z,m_He)/8) + C_D_caja*L_x*L_y*rho_atm(z-R_globo(z,m_He)-l_caja_globo));

z_int = z_0:dz:z_f;
int1 = regla_de_simpson(dz, z_int, @(z)E(z,m_He)-F_g_caja_z(z));
int2 = regla_de_simpson(dz, z_int, F_roz_norm_z);

M = M_caja+m_He+M_paraca1+M_paraca2+M_globo;    % Masa total durante el ascenso (kg)

% Cálculo de la aceleración promedio (m/s^2)
a = (int1 - G*M_T*(M_globo+m_He)*(1/(R_T+z_0)-1/(R_T+z_f)))/(M*(z_f-z_0) - int2);

path(ruta)

% Cálculo de las velocidades medias (en m/s) para el MRUA, promediando 
% respecto a la altitud y respecto al tiempo.
% Se hace uso de la propiedad: <v>_z(t) = <v*dz_dt>_t/<dz_dt>_t
v_prom_z = 2 * sqrt(2 * a * (z_f-z_0)) / 3;
v_prom_t = sqrt(a * (z_f-z_0) / 2);


% NOTA: para los vuelos que se han simulado, la media aritmética de los
% valores de la velocidad simulada está siempre entre v_prom_z y v_prom_t,
% aunque normalmente se queda más cerca del valor v_prom_t.