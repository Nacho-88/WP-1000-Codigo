function M_He = calcularMasaHelio(R_exp, R_0, z_exp, T_amb, P_amb, validacion)
load constantes.mat R R_prima

if validacion
    % Obtener la masa a partir de unas condiciones iniciales
    volumen = (4/3) * pi * (R_0^3);
    densidad = P_amb / (R_prima*(T_amb+273.15));
else
    % Obtener la masa necesaria para alcanzar una cierta altura
    volumen = (4/3) * pi * (R_exp^3);
    densidad = (R/R_prima) * rho_atm(z_exp);
end

M_He = volumen * densidad;

end

