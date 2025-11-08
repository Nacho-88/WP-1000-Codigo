function [P_n] = generar_P_n()

load constantes.mat z_n R_T P_0

N=length(z_n);

P_n=zeros(1,N);

P_n(1)=P_0;

for i=2:N

    P_n(i)= P(z_n(i-1)*(R_T/(R_T-z_n(i-1))));

end


end