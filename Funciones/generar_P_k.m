function [P_k] = generar_P_k()

load constantes.mat z_star P_n

n=n_capa(z_star);

N=lentgh(P_n);

P_k=zeros(1,N);

P_0=P(z_star);

P_k(1)=P_0;

for i=2:N

  P_k(i)= P_n(i+n);

end

end