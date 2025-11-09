% F_g_caja.m
% ---------------------------------------------------------------
% Cálculo de la fuerza gravitatoria del payload (caja)
% 
% Fórmula:
%   F_g(z) = (G * M_T * m_caja) / (R_T + z - r_globo(z) + l_c - L_z/2)^2
%
% Ámbito: válido tanto en ascenso como en descenso
%
% Inputs:
%   z         -> altitud real [m]
%   dz_dt     -> velocidad vertical [m/s] (actualmente no se usa)
%   isFalling -> booleano (true si está descendiendo, false si ascendiendo)
%
% Output:
%   F         -> fuerza gravitatoria [N]
%
% Dependencias (definidas en el main):
%   G, M_T, R_T, m_caja, l_c, L_z
%   r_globo(z)
% ================================================================

function F = F_g_caja(z, m_He, isFalling) 

load constantes.mat G M_T R_T M_caja L_z l_caja_globo
    
    % Distinción entre ascenso y descenso con el nuevo parámetro
    if isFalling
        %DESCENSO
       
        F = G *M_T * M_caja / ((R_T + z - L_z/2)*(R_T + z + L_z/2));
    else
        %ASCENSO
       
        F = G *M_T * M_caja / ((R_T + z - (R_globo(z, m_He) + l_caja_globo + L_z))*(R_T + z - (R_globo(z, m_He) + l_caja_globo));
    end

end




