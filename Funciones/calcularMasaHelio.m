function M_He = calcularMasaHelio(R_exp, R_0, z_exp, validacion)
% % 1. Término del Volumen (4/3 * pi * R^3)
% termino_volumen = (4/3) * pi * (R_exp^3);
% 
% % 2. Término del ratio de presión elevado al exponente
% n_exp = n_capa(z_exp);
% 
% load constantes.mat P_n T_k z_star R_prima
% Pm_zexp = P(z_exp);
% P_np = P_n(n_exp);
% termino_presion = (Pm_zexp / P_np)^(3/5);
% 
% % 3. Término del factor final
% n_star = n_capa(z_star);
% T_ref = T_k(n_exp-n_star);
% termino_factor = (P_np / (R_prima*T_ref));
% 
% % --- Cálculo final ---
% % M_He = (Término 1) * (Término 2) * (Término 3)
% M_He = termino_volumen * termino_presion * termino_factor;

load constantes.mat R R_prima L_z l_caja_globo

if validacion
    % Obtener la masa a partir de unas condiciones iniciales
    volumen = (4/3) * pi * (R_0^3);
    densidad = (R/R_prima) * rho_atm(1026 + L_z + l_caja_globo + R_0);
else
    % Obtener la masa necesaria para alcanzar una cierta altura
    volumen = (4/3) * pi * (R_exp^3);
    densidad = (R/R_prima) * rho_atm(z_exp);
end

M_He = volumen * densidad;

end

