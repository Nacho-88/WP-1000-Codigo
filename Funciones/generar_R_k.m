function [R_k] = generar_R_k(m_He)

load constantes.mat z_star P_k 

N=20;

R_k=zeros(1,N);

R_0=R_globo(z_star, m_He); 

R_k(1)=R_0;

for i=2:N

    R_k(i) = R_k(i-1)*(P_k(i-1)/ P_k(i))^(1/5);

end

end