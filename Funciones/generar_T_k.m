function [T_k] = generar_T_k()

load constantes.mat z_star P_k

N=length(P_k);

T_0=T_atm(z_star); 

T_k=zeros(1,N);

T_k(1)=T_0;

for i=2:N

    T_k(i)= T_k(i-1)*(P_k(i)/P_k(i-1))^(2/5);


end

end