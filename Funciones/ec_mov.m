function [dw_dt] = ec_mov(t, w, m_He)
% EC_MOV: función que define el sistema de ecuaciones diferenciales: f(y,t) = dy/dt
%
% -'t' es el tiempo (en segundos).
%
% -'w' es un vector (columna) definido como: [dz_dt, z]' donde z es la
% altitud real (en metros) del globo/payload y dz_dt es su derivada temporal (m/s).
%
% -'m_He' es la masa de Helio que contiene el globo (en kilogramos).
%
%
% -'dw_dt' es la derivada temporal de 'w', que en este caso es igual a [d^2z_dt^2, dz_dt]'.


%% Comprobación valores de inputs
if m_He<=0
    warning("Error en el parámetro de la función 'ec_mov': la masa debe ser positiva.")
end

if w(2)<0
    warning("Error en el parámetro de la función 'ec_mov': la altitud (w(2)) debe ser no negativa.")
end

if t<0
    warning("Error en el parámetro de la función 'ec_mov': el tiempo debe ser no negativo.")
end


%% Cargar constantes necesarias
% load ../constantes.mat m_latex m_caja m_paraca_1 m_paraca_2 % Versión para ejecutar desde la carpeta "Funciones"
load constantes.mat m_latex m_caja m_paraca_1 m_paraca_2 % Versión para ejecutar desde el main.


%% Definición de la ecuación del movimiento

% Comprobamos si estamos en ascenso o descenso (dz_dt>0 ó dz_dt<0 resp.)

if w(1)>0 | t==0    % Cuando t=0 entonces dz_dt=0, en cualquier otro caso donde dz_dt=0 no estamos en ascenso

    dw_dt = [(E(w(2), m_he) + F_roz(w(2),w(1)) - F_g(w(2), w(1), m_He))/(m_He + m_latex + m_caja + m_paraca_1 + m_paraca_2), ...
            w(1)]';

elseif w(1)<=0 & t>0    % Si dz_dt=0 necesitamos que t>0 para estar en el 
                        % punto final o en la máxima altitud (pasamos de 
                        % ascenso a descenso, por lo que para calcular los 
                        % siguientes puntos usamos la ecuación de descenso)

    dw_dt = [(F_roz(w(2), w(1)) - F_g(w(2), w(1), m_He))/(m_caja + m_paraca_1 + m_paraca_2), ...
            w(1)]';

else    % No deberíamos entrar aquí

    warning("Situación sin sentido. (ec_mov)")

    dw_dt = NaN;

end


end

