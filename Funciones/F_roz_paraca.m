function F_paraca = F_roz_paraca(z, dz_dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la fuerza de rozamiento debido a ambos paracaídas a partir de la altitud (real) y la velocidad de descenso
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg), dz_dt (velocidad ascenso/descenso en m/s)
% Salida: Fuerza: F (medida en Newtons)
% Constante: C_D_paraca_1, C_D_paraca_2, l_c (m), D_1 (diámetro paracaídas 1, medido en metros), D_2 (m)
% Depende del Radio del globo (R_globo), Densidad atmosfera (rho_Atm) y Presión atmósfera (P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coeficientes de arrastre
C_D_paraca_1 = 0.97;
C_D_paraca_2 = 1.17;

% Areas paracaidas
A1 = 3.25;
A2 = 1.13;

% Longitud de cuerdas
l_caja_paraca1  = 12.0;    % Caja – Paracaídas 1
l_paraca1_paraca2 = 15.0;  % Paracaídas 1 – Paracaídas 2
% l_paraca2_globo = 10.0;       % Paracaídas 2 – Globo

% Densidad de la atmosfera
z_paraca1 = z + l_caja_paraca1;
rho1 = rho_atm(z_paraca1);

z_paraca2 = z + l_caja_paraca1 + l_paraca1_paraca2;
rho2 = rho_atm(z_paraca2);

% Solo descenso
if dz_dt >= 0
    F_paraca = 0;
    return;
end

% Fuerza de rozamiento
F_paraca = (0.5 * (dz_dt.^2) * rho1 * C_D_paraca_1*A1) + (0.5 * (dz_dt.^2) * rho2 *C_D_paraca_2*A2);



end

