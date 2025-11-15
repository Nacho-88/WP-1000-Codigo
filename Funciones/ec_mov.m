function [dw_dt] = ec_mov(t, w, m_He, R_exp)
% EC_MOV: función que define el sistema de ecuaciones diferenciales: f(t,y) = dy/dt
%
% -'t' es el tiempo (en segundos).
%
% -'w' es un vector (columna) definido como: [dz_dt, z]' donde z es la
% altitud real (en metros) del globo/payload y dz_dt es su derivada temporal (m/s).
%
% -'m_He' es la masa de Helio que contiene el globo (en kilogramos).
%
% -'R_exp' es el radio con el que el globo explota (en metros).
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

if R_exp<=0
    warning("Error en el parámetro de la función 'ec_mov': el radio de explosión debe ser positivo.")
end


%% Cargar constantes necesarias
% load ../constantes.mat m_globo m_caja m_paraca_1 m_paraca_2 l_caja_globo L_z % Versión para ejecutar desde la carpeta "Funciones"
load constantes.mat M_globo M_caja M_paraca1 M_paraca2 l_caja_globo L_z % Versión para ejecutar desde el main.


%% Definición de la ecuación del movimiento

% Comprobamos si estamos en ascenso o descenso (dz_dt>0 ó dz_dt<0 resp.)

if (w(1)>0 | t==0) & (R_globo(w(2),m_He)<R_exp)  % Cuando t=0 entonces dz_dt=0, en cualquier otro caso donde dz_dt<=0 (o el radio no sea inferior al de explosión) no estamos en ascenso
    
    isFalling = false;  % Estamos en ascenso.
    dw_dt = [(E(w(2), m_He) + F_roz(w(2),w(1),m_He,isFalling) - F_g(w(2),m_He,isFalling))/(m_He + M_globo + M_caja + M_paraca1 + M_paraca2), ...
            w(1)]';

elseif (w(1)<=0 & t>0)  % Si dz_dt=0 necesitamos que t>0 para estar en el 
                        % punto final o en la máxima altitud (pasamos de 
                        % ascenso a descenso, por lo que para calcular los 
                        % siguientes puntos usamos la ecuación de descenso)

    isFalling = true;   % Estamos en descenso.
    dw_dt = [(F_roz(w(2), w(1), m_He, isFalling) - F_g(w(2), m_He, isFalling))/(M_caja + M_paraca1 + M_paraca2), ...
            w(1)]';

elseif (w(1)>0) & (R_globo(w(2),m_He)>=R_exp)    % Momento de explosión del globo

    % Cambiamos la referencia de la altitud del centro del globo al centro
    % del payload:
    w(2) = w(2) - R_globo(w(2),m_He) - l_caja_globo - L_z/2;

    % Guardamos instante y altitud en el que detectamos la explosión. Como
    % podemos estar en varios instantes de tiempo en la transición,
    % guardamos todos los tiempos en los que nos encontramos en esta
    % situación, de tal forma que el primero corresponde con el tiempo de
    % explosión.
    load parametros.mat t_exp z_exp
    
    if ~exist("t_exp","var")
        t_exp = t;
        save parametros.mat t_exp -append
    end

    if ~exist("z_exp","var")
        z_exp = w(2);
        save parametros.mat z_exp -append
    end
    

    isFalling = true;   % Pasamos de ascenso a descenso.
    dw_dt = [(F_roz(w(2), w(1), m_He, isFalling) - F_g(w(2), m_He, isFalling))/(M_caja + M_paraca1 + M_paraca2), ...
            w(1)]';

else    % No deberíamos entrar aquí

    warning("Situación sin sentido. (ec_mov)")

    dw_dt = [NaN, NaN]';

end


end

