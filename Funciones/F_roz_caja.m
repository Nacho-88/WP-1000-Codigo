function F_caja = F_roz_caja(z, dz_dt, m_he)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la fuerza de rozamiento del payload a partir de la altitud (real) y la velocidad de ascenso o descenso
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg), dz_dt (velocidad ascenso/descenso en m/s)
% Salida: Fuerza: F (medida en Newtons)
% Constante: C_D_semiesfera, A_caja (Lx*Ly)
% Depende del Radio del globo (R_globo), Densidad atmosfera (rho_Atm) y Presión atmósfera (P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constantes
C_D_caja = 1.57439;     % Coeficiente de arrastre caja            
l_caja_globo = 37;      % Longitud (m) de la cuerda entre globo y payload

% Área caja
A_caja = 0.098;     % [m^2]

% Radio del globo
R_g = R_globo(z, m_he);

% Densidad de la atmosfera
z_caja = z-R_g-l_caja_globo;
rho = rho_atm(z_caja);

% Signo de la velocidad: ascenso (+1), descenso (-1)
if dz_dt > 0
    signo_velo = 1;
elseif dz_dt < 0
    signo_velo = -1;
end

% Fuerza de rozamiento
F_caja = 0.5 * (dz_dt.^2) * C_D_caja * A_caja * rho * (-signo_velo);


end

