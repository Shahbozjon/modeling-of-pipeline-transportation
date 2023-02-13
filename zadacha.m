%   ishodniye danniye
p1 = 6*10^6;
rho = 860;
visc = 12*10^(-6);
delta = 0.1/1000;
d = 0.4234;
Q = 800/3600;
L = 100*10^3;
z = [41, 14, 53, 63, 354, 24, 64, 12, 46, 13, 52];
%--------------
x = 0:10000:L;
v = 4*Q/pi()/d^2;
h = lambda(v,d,visc,delta)*L/d*v^2/2/9.81;
H1 = p1/rho/9.81 + z(1);
H2 = H1 - h;

%postroenie grafikov
plot([x(1),x(11)],[H1,H2],'LineWidth',2)
hold on
plot(x,z,'r','LineWidth',2)
hold off