function M_He = calcularMasaHelio(R_exp, z_exp)
% 1. Término del Volumen (4/3 * pi * R^3)
termino_volumen = (4/3) * pi * (R_exp^3);

% 2. Término del ratio de presión elevado al exponente
n = n_capa(z_exp);
[P_n] = generar_P_n();
[P_k] = generar_P_k();
Pm_zexp = P_k(n);
P_np = P_n(n);
termino_presion = (Pm_zexp / P_np)^(3/5);

% 3. Término del factor final
[T_k] = generar_T_k() ;
T_ref = T_k(n);
termino_factor = (P_np / (R_prime*T_ref));

% --- Cálculo final ---
% M_He = (Término 1) * (Término 2) * (Término 3)
M_He = termino_volumen * termino_presion * termino_factor;

end
