%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOMBRES VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%altura inicial z0
%altura z
%gravedad g
%Empuje E
%masa helio m_He
%temperatura atmosferica T_atm
%temperatura globo T_globo
%presion atmosferica P_atm (Pa)
%presion globo P_globo
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
r_globo=1;
M_T=5.97e24;
G=6.67e-11;
g=9.81;
P_atm=101325;
R_T=6.37e6;
c_d=0.75;
m_He=2;
dpar=1.64;
m_atomicaHe=4;
Cd=0.3;
v=5.33;

%Fórmula temperatura atmosférica
Tp = @(z) Tn + (z - zn) * (Tnp1 - Tn) / (znp1 - zn);

altitud_m = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40] * 1000;
temperatura_K = [15, 2, -10.98, -23.96, -36.93, -49.9, -56.5, -56.5, -56.5, -56.5, ...
                 -56.5, -54.5, -52.5, -50.5, -48.62, -46.64, -44.66, -39.41, -33.87, -28.33, -22.8] + 273.15;


