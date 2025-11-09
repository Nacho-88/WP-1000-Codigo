function F_caja = F_roz_caja(z, dz_dt, m_He,isFalling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la fuerza de rozamiento del payload a partir de la altitud (real) y la velocidad de ascenso o descenso
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg), dz_dt (velocidad ascenso/descenso en m/s)
% Salida: Fuerza: F (medida en Newtons)
% Constante: C_D_semiesfera, A_caja (Lx*Ly)
% Depende del Radio del globo (R_globo), Densidad atmosfera (rho_Atm) y Presión atmósfera (P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Llamar constantes
load constantes C_D_caja l_caja_globo L_z A_caja

% Radio del globo
R_g = R_globo(z, m_He);

% Fuerza de rozamiento: ascenso, descenso
if ~isFalling
    z_caja = z-R_g-l_caja_globo;
    signo_velo = +1;
else
    z_caja = z+(L_z/2);
    signo_velo = -1;
end

rho = rho_atm(z_caja);

% Fuerza de rozamiento
F_caja = 0.5 * (dz_dt.^2) * C_D_caja * A_caja * rho * (-signo_velo);


end

