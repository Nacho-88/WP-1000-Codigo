function [P_k] = generar_P_k()

load constantes.mat z_star P_n

n=n_capa(z_star);

N=length(P_n);

P_k=zeros(1,N);

P_0=P(z_star);

P_k(1)=P_0;

for i=2:N

  P_k(i)= P_n(i+n-1);

  if (i+n-1)==N   % El índice i+n-1 no debe superar la longitud del array
      return
  end

end

end