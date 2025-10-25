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
R_T = 6.37e6;   % m
g_0 = 9.81;     % m / s^2
P_0 = 101325;   % Pa

T_n = [15, -56.5, -56.5, -44.5, -2.5, -2.5, -58.5]; % ºC
T_n = T_n + 273.15.*ones(size(T_n));        % K
z_n = 1e3 .* [0, 11, 20, 32, 47, 51, 71];   % m


r_globo=1;
c_d=0.75;
m_He=2;
dpar=1.64;
m_atomicaHe=4;
Cd=0.3;
v=5.33;

% %Fórmula temperatura atmosférica
% Tp = @(z) Tn + (z - zn) * (Tnp1 - Tn) / (znp1 - zn);
% 
% altitud_m = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40] * 1000;
% temperatura_K = [15, 2, -10.98, -23.96, -36.93, -49.9, -56.5, -56.5, -56.5, -56.5, ...
%                  -56.5, -54.5, -52.5, -50.5, -48.62, -46.64, -44.66, -39.41, -33.87, -28.33, -22.8] + 273.15;

save constantes G R_T M_T g_0 P_0 T_n z_n



%% Resolución ecuación del movimiento



%% Gráficas y datos importantes a comparar
