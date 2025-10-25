function [n] = n_capa(z)
% Función que, dada una altitud real, devuelve el índice n correspondiente 
% a la capa de la atmósfera del modelo.
%
% -'z' es la altitud real medida en metros.
%
%
% -'n' es el índice que indica en qué capa del modelo de la atmósfera se
%  encuentra la altitud dada.



% % Comprobación valores de inputs:
% if z<0
%     error("Error en el parámetro de la función 'n_capa': la altitud debe ser positiva.")
% end


% Cargar valores de constantes necesarias
load ../constantes.mat T_n z_n R_T


% Inicio del código funcional
n = 0;

for i=1:length(z_n)

    if z<z_n(i+1)*R_T/(R_T-z_n(i+1)) & z>=z_n(i)*R_T/(R_T-z_n(i))
        n = i-1;  % Encontramos la capa buscada. (El índice 'i' está desplazado una unidad por encima porque MATLAB empieza en 1 en vez de en 0)
        return  % Hemos terminado.
    end

end

% Si llegamos a esta parte del código es porque no se ha encontrado la capa
% (estamos a una altitud superior)
n = length(z_n);

end

