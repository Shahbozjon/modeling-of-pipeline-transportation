%Neustanovivsheesa techeniye nefti
clc; clear
t = 0;
p_nach = 6*10^6; p_kon = 2*10^6;
rho = 860;
visc = 12*10^(-6);
delta = 0.1/1000;
d = 1;
Q = -4800/3600;
L = 100*10^3;

c = 1000; %skorost' zvuka
N = 100; dx = L/N; dt = dx/c;
x = 0:dx:L;

v0 = zeros(1,N+1); v0(:) = 4*Q/pi()/d^2;
p0 = zeros(1,N+1);

for i = 1:N+1
   p0(i) =  p_nach - (p_nach-p_kon)/ L * dx*(i-1);
end

for j = 1:1000
t = t + dt;
if t > 100
    p_kon = 5*10^6;
end
if t > 120
    p_nach = 3*10^6;
end
p1 = zeros(1,N+1); v1 = zeros(1,N+1);
%levaya granica
i = 1;
p1(i) = p_nach;
InvB = p0(i+1) - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
v1(i) = (p1(i) - InvB) / (rho*c);

%tsentralnaya chast truby
for i = 2:N
   InvA = p0(i-1) + rho*c*v0(i-1) - lambda(v0(i-1),d,visc,delta)*rho*v0(i-1)*abs(v0(i-1))/2/d*c*dt;
   InvB = p0(i+1) - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
   p1(i) = (InvA+InvB)/2;
   v1(i) = (InvA-InvB)/(2*rho*c);
end

%pravaya granica
i = N+1;
p1(i) = p_kon;
InvA = p0(i-1) + rho*c*v0(i-1) - lambda(v0(i-1),d,visc,delta)*rho*v0(i-1)*abs(v0(i-1))/2/d*c*dt;
v1(i) = (InvA - p1(i)) / (rho*c);

p0 = p1; v0 = v1;

subplot(2,1,1)
    plot(x,p1,'LineWidth',2);
    grid on
    axis([0 L 0 80*10^5])
    xlabel('length (km)','fontsize',16)
    ylabel('Presssure (Pa)','fontsize',16)
    title(sprintf('time = %1.1f',t))
     
    subplot(2,1,2)
    plot(x,v1,'LineWidth',2)
    grid on
    axis([0 L -6 6])
    xlabel('length (km)','fontsize',16)
    ylabel('Velocity (m/s)','fontsize',16)
    shg
    pause(0.0001)
end