function dx = Fitzhugh_Nagumo(t,x,a,b,I_ext,tau)

dx =[x(1) - x(1)^3/3 - x(2) + I_ext;
     1/tau*(x(1) + a - b*x(2))];