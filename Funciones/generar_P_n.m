function [P_n] = generar_P_n()

load constantes.mat z_n P_0 T_n R g_0

N=length(z_n);

P_n=zeros(1,N);

P_n(1)=P_0;

for i=2:N

    if T_n(i)~=T_n(i-1)
        base = 1 + (T_n(i)-T_n(i-1))/(T_n(i-1));
        exponente = -g_0*(z_n(i)-z_n(i-1))/(R*(T_n(i)-T_n(i-1)));
    else
        base = exp(1);  % Número de Euler
        exponente = -g_0*(z_n(i)-z_n(i-1))/(R*T_n(i-1));
    end

    P_n(i) = P_n(i-1)*base^exponente;

end


end