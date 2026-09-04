function dYdt = Jansen_Rit(t,Y,parameters)

dYdt = zeros(6,1);

%parameters setting
C      = parameters(1);
A      = parameters(2);
B      = parameters(3);
a      = parameters(4);
b      = parameters(5);
V_0    = parameters(6); 
nu_max = parameters(7);
r      = parameters(8); 
p      = parameters(9);

% Jansen-Rit model equations
y_1 = Y(1); y_2 = Y(2); y_3 = Y(3); y_4 = Y(4); y_5 = Y(5); y_6 = Y(6);
dYdt(1) = y_4;
dYdt(2) = y_5;
dYdt(3) = y_6;
dYdt(4) = A*a*(nu_max/(1+exp(r*(V_0-(y_2-y_3))))) - 2*a*y_4 - a^2*y_1;
dYdt(5) = A*a*(p + 0.8*C*(nu_max/(1+exp(r*(V_0-(C*y_1)))))) - 2*a*y_5 - a^2*y_2;
dYdt(6) = B*b*0.25*C*(nu_max/(1+exp(r*(V_0-(0.25*C*y_1))))) - 2*b*y_6 - b^2*y_3;

end