function [P_b, delta_P, Nozzle_Lift] = llenado(m_He, T_amb)

load constantes.mat R_prima
P_b0 = 200;         % Bar
V_b = 50*10^-3;     % m^3

tz = T_amb + 273.15;

delta_P = (m_He*R_prima*tz/V_b) * 10^(-5);  % bar
P_b = P_b0 - delta_P;                       % bar

Nozzle_Lift = E((1026), m_He) - F_g_globo((1026), m_He);    % N

end