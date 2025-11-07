function F_total = F_roz(z, dz_dt, m_he)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcula la fuerza de rozamiento con el aire (tanto en ascenso como en descenso) a partir de la altitud (real/geométrica). Tanto ascenso como descenso
% Entrada: z (altitud geometrica), m_he (amsa de Helio en kg), dz_dt (velocidad ascenso/descenso en m/s)
% Salida: Fuerza: F (medida en Newtons)
% Depende de la Fuerza rozamiento globo, Fuerza rozamiento caja, Fuerza rozamiento paracaídas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dz_dt > 0                                   % Ascenso
    
    F_globo = F_roz_globo(z, dz_dt, m_he);
    F_caja  = F_roz_caja(z, dz_dt, m_he);
    F_total = F_globo + F_caja;

else                                           % Descenso
   
    F_paraca = F_roz_paraca(z, dz_dt);
    F_caja   = F_roz_caja(z, dz_dt, m_he);
    F_total  = F_paraca + F_caja;
    
end


end

