%% Descripción
% Queremos hacer una estimación de la altitud y velocidad de ascenso a
% partir de promedios temporales. Para ello aproximamos que la aceleración
% de ascenso en un cierto intervalo de tiempo se mantiene constante y, por
% tanto, es igual a su promedio temporal. A partir de este valor de la
% aceleración se deduce la velocidad promedio hasta ese instante de tiempo.


%% Cálculo de la aceleración

load parametros.mat t_exp
z_0 = 1026;     % m
dz = 10;        % m
dt2 = 6480;     % s

R_exp = 12.4/2; % m
R_0 = 2.25/2;   % m
z_exp = 37550;  % m
validacion = true;

t_f = 1800;     % s

ruta = addpath('Funciones');
m_He = calcularMasaHelio(R_exp, R_0, z_exp, validacion);
path(ruta)

%% Cálculo directo del promedio de la aceleración hasta explosión
%z_f = z_0 + 300;
z_f = z_exp;

ruta = addpath('Funciones');
load constantes.mat G M_T R_T M_caja l_caja_globo L_z C_D_globo C_D_caja L_x L_y M_paraca1 M_paraca2 M_globo
F_g_caja_z = @(z) G*M_T*M_caja/((R_T+z-R_globo(z,m_He)-l_caja_globo-L_z)*(R_T+z-R_globo(z,m_He)));
F_roz_norm_z = @(z) -(z-z_0)*(pi*C_D_globo*(R_globo(z,m_He))^2*rho_atm(z+3*R_globo(z,m_He)/8) + C_D_caja*L_x*L_y*rho_atm(z-R_globo(z,m_He)-l_caja_globo));

z_int = z_0:dz:z_f;
int1 = regla_de_simpson(dz, z_int, @(z)E(z,m_He)-F_g_caja_z(z));
int2 = regla_de_simpson(dz, z_int, F_roz_norm_z);

M = M_caja+m_He+M_paraca1+M_paraca2+M_globo;    % kg

a = (int1 - G*M_T*(M_globo+m_He)*(1/(R_T+z_0)-1/(R_T+z_f)))/(M*(z_f-z_0) - int2);

path(ruta)

v_prom_exp_z = 2 * sqrt(2 * a * (z_f-z_0)) / 3;
v_prom_exp_t = sqrt(a * (z_f-z_0) / 2);


%% Cálculo aceleración promedio hasta t=t_f
ec_t = @(a) ecuacion_t(a, t_f, z_0, dt2, m_He);


% Llamar bisección para resolver (para 'a') ec_t(a) = 0 
a_t = Met_Biseccion(ec_t, 1e-8, 20, 0, 1e-5);

% Calcular la velocidad promedio a partir del valor de 'a':
% Como estamos promediando respecto a la altitud en vez del tiempo, hacemos
% uso de la siguiente identidad: <F>_x(t) = <F*dx/dt>_t/<dx/dt>_t donde en
% este caso F=v=dz/dt y x(t)=v(t)=a*t^2/2
% Por lo que nos queda: <v>_z = <v^2>_t/<v>_t = 2*a*t_f/3
v_prom_z = 2 * a_t * t_f/3;
v_prom_t = a_t * t_f/2;


%% Cálculo aceleración promedio hasta explosión (z=z_exp)
ec_z = @(a) ecuacion_z(a, z_exp, z_0, dz, m_He);

% Llamar bisección para resolver (para 'a') ec_z(a) = 0
a_z = Met_Biseccion(ec_z, 1e-8, 1e-3, 0, 1e-5);

% Calcular la velocidad promedio a partir del valor de 'a': <v>_z =
% 2*a*t_exp/3  =>  <v>_z = 2*sqrt(2*(z_exp-z_0)*a)/3
v_prom_exp_z = 2 * sqrt(2 * a_z * (z_exp-z_0)) / 3;
v_prom_exp_t = sqrt(a_z * (z_exp-z_0) / 2);


%% Funciones auxiliares
function output = ecuacion_z(a, z_f, z_0, dz, m_He)

    load constantes.mat M_T M_caja M_paraca1 M_paraca2 G R_T M_globo

    ruta = addpath('Funciones');


    x_inf = z_0;
    x_sup = z_f;
    dx = dz;

    x_int = x_inf:dx:x_sup;

    E_z2 = @(z) E_t(z_0, a, (z-z_0)/a, m_He);
    F_roz_z2 = @(z) F_roz_t(z_0, a, (z-z_0)/a, m_He);
    F_g_caja_z2 = @(z) F_g_caja_t(z_0, a, (z-z_0)/a, m_He);

    % E_int = regla_de_simpson(dx, x_int, E_z2);
    % F_roz_int = regla_de_simpson(dx, x_int, F_roz_z2);
    % F_g_int = G*M_T*(M_globo+m_He)*(1/(R_T+z_0)-1/(R_T+z_f)) + regla_de_simpson(dx, x_int, F_g_caja_z2);

    int_prom = regla_de_simpson(dx, x_int, @(z) E_z2(z)+F_roz_z2(z)-F_g_caja_z2(z));

    path(ruta);
    
    % output = E_int + F_roz_int - F_g_int - (M_caja+m_He+M_paraca1+M_paraca2)*(a*(z_f-z_0));
    output = int_prom - G*M_T*(M_globo+m_He)*(1/(R_T+z_0)-1/(R_T+z_f)) - (M_caja+m_He+M_paraca1+M_paraca2)*a*(z_f-z_0);

end


function output = ecuacion_t(a, t_f, z_0, dt2, m_He)

    load constantes.mat G M_T M_globo R_T M_caja M_paraca1 M_paraca2

    ruta = addpath('Funciones');


    x_inf = 0;
    x_sup = t_f^2/2;
    dx = dt2;

    x_int = x_inf:dx:x_sup;

    E_t2 = @(t) E_t(z_0, a, t, m_He);
    F_roz_t2 = @(t) F_roz_t(z_0, a, t, m_He);
    F_g_caja_t2 = @(t) F_g_caja_t(z_0, a, t, m_He);

    E_int = regla_de_simpson(dx, x_int, E_t2);
    F_roz_int = regla_de_simpson(dx, x_int, F_roz_t2);
    F_g_int = G*M_T*(M_globo+m_He)*t_f^2/(2*(R_T+z_0)*(R_T+z_0+a*t_f^2/2)) + regla_de_simpson(dx, x_int, F_g_caja_t2);


    path(ruta);
    
    output = E_int + F_roz_int - F_g_int - (M_caja+m_He+M_paraca1+M_paraca2)*(a*t_f^2/2);

end


function Empuje = E_t(z_0, a, t, m_He)

    Empuje = E(z_0+a*t, m_He);

end


function Fuerza = F_roz_t(z_0, a, t, m_He)

    load constantes.mat C_D_caja C_D_globo l_caja_globo L_x L_y

    R_g = R_globo(z_0+a*t, m_He);
    Fuerza = -a^2*t/2*(pi*(R_g)^2*C_D_globo*rho_atm(z_0+a*t+3*R_g/8) + C_D_caja*L_x*L_y*rho_atm(z_0+a*t-R_g-l_caja_globo));

end


function Fuerza = F_g_caja_t(z_0, a, t, m_He)

    load constantes.mat l_caja_globo G M_T L_z M_caja R_T

    z_aux = R_T + z_0 + a*t - R_globo(z_0+a*t, m_He) - l_caja_globo;

    Fuerza = G*M_T*M_caja/(z_aux*(z_aux-L_z));

end