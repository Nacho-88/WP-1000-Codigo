function [outputArg1,outputArg2] = untitled(z0,T_atm,P_atm,g,G,R_T,M_T,m_paracaidas,c_d,dpar,m_He,altitud_m,temperatura_K,m_atomicaHe,r_globo,E,M,Cd)
for i=1: 1000000
    %Calculo presión atmosférica
if temperatura_K(i) == temperatura_K(i+1)
    P_atm = @(z) P_atm * exp(((-g * z - altitud_m(i)) / (R_T * temperatura_K(i))));
else
    P_atm = @(z) P_atm * (1 + ((z - altitud_m(i)) * (temperatura_K(i+1) - temperatura_K(i))) / ((altitud_m(i+1) - altitud_m(i)) * temperatura_K(i))) ^ ...
        ((-g * (altitud_m(i+1) - altitud_m(i))) / (R_T * (temperatura_K(i+1) - temperatura_K(i))));
end

F_roz = 1/2*Cd*pi*(2r_globo(z)^2)/4*rho
if z<12000
    r_globo = @(z) (3*m_He*(m_atomicaHe/0.0289644)*r_globo(i-1)*temperatura_K)/(4*pi*P_atm(z))^1/3;
else z >= 12000
    r_globo = @(z) Rk*(Pk/P_atm(z))^1/5
end
E = @(z)(P_atm*(z+sqrt((r_globo^2)-r^2))-P_atm*(z-sqrt((r_globo^2)-r^2))*r*v

dx=0.5;

integral = regla_de_simpson(dx, x, f)


dt=0.5;
t0=0;
tf=100000;
y0=1;

h= @(z) E/M + Froz/M- (G*M_T)/(R_T+z)^2

[times, y_aprox] = Metodo_RK4(dt, t0, tf, f, y0)

if(P_atm>P_globo)

    return;

end


end