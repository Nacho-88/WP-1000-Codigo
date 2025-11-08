
% ================================================================
% F_g.m
% ---------------------------------------------------------------
% Fuerza gravitatoria TOTAL sobre la unidad de vuelo:
%
%   F_g(z) = F_g_caja(z, dz_dt, isFalling) + F_g_globo(z, dz_dt, m_He)
%            + G*M_T*m_paraca_1 / (R_T + z - r_globo(z) - l_1)^2
%            + G*M_T*m_paraca_2 / (R_T + z - r_globo(z) - l_2)^2
%
% Prototipo:
%   F = F_g(z, dz_dt, m_He, isFalling)
%
% Inputs:
%   z         : altitud [m]
%   dz_dt     : velocidad vertical [m/s]
%   m_He      : masa de Helio [kg]
%   isFalling : booleano -> true si está en descenso, false si en ascenso
%
% Output:
%   F         : fuerza gravitatoria total [N]
%
% Dependencias (definidas en el main):
%   Constantes: G, M_T, R_T, m_paraca_1, m_paraca_2, l_1, l_2
%   Funciones : F_g_caja(z, dz_dt, isFalling), F_g_globo(z, dz_dt, m_He), r_globo(z)
% ================================================================

function F = F_g(z, dz_dt, m_He, isFalling)

    %Caja: distingue ascenso/descenso
    F_caja  = F_g_caja(z, dz_dt, isFalling);

    %Globo: SIEMPRE se incluye (no distingue ascenso/descenso)
    F_globo = F_g_globo(z, dz_dt, m_He);

    %Fuerza de los paracaídas
    F_paraca1 = (G * M_T * m_paraca_1) ./ (R_T + z - r_globo(z) - l_1).^2;
    F_paraca2 = (G * M_T * m_paraca_2) ./ (R_T + z - r_globo(z) - l_2).^2;

    %Fuerza total
    F = F_caja + F_globo + F_paraca1 + F_paraca2;
end

