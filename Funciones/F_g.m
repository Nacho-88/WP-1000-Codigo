
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

function F = F_g(z, m_He, isFalling)

load constantes.mat G M_T R_T M_paraca1 M_paraca2 l_caja_paraca1 l_paraca1_paraca2 l_paraca2_globo L_z

    %Caja: distingue ascenso/descenso
    F_caja  = F_g_caja(z, m_He, isFalling);

    %Globo: SIEMPRE se incluye (no distingue ascenso/descenso)
    F_globo = F_g_globo(z, m_He);

    %Fuerza total
    if isFalling
        %Fuerza de los paracaídas en descenso
        F_paraca1 = (G * M_T * M_paraca1) ./ (R_T + z + L_z/2 + l_caja_paraca1).^2;
        F_paraca2 = (G * M_T * M_paraca2) ./ (R_T + z + L_z/2 + l_caja_paraca1 + l_paraca1_paraca2).^2;

        F = F_caja + F_paraca1 + F_paraca2;
    else
        %Fuerza de los paracaídas en ascenso
        F_paraca1 = (G * M_T * M_paraca1) ./ (R_T + z - R_globo(z,m_He) - (l_paraca1_paraca2 + l_paraca2_globo)).^2;
        F_paraca2 = (G * M_T * M_paraca2) ./ (R_T + z - R_globo(z,m_He) - l_paraca2_globo).^2;

        F = F_caja + F_globo + F_paraca1 + F_paraca2;
    end
end

