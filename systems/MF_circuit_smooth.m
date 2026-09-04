function dx = MF_circuit_smooth(t,x,L_1,L_2,C_1,R,V_S,V_th,C_2,beta,k)

Gamma = @(y) y*(1+tanh(k*y))/2;

dx1 = (V_S-x(1))/(R*C_1) - (x(3)+x(4))/C_1;
dx2 = (x(4)-beta*Gamma(x(3))*tanh(x(2)/(2*V_th)))/C_2;
dx3 = (x(1)-V_th)/L_1;
dx4 = (x(1)-x(2))/L_2;

dx =[dx1;dx2;dx3;dx4];
