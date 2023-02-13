%Neustanovivsheesa techeniye nefti
clc; clear
t = 0;
p_nach = 6*10^6; p_kon = 2*10^6;
rho = 860;
visc = 12*10^(-6);
delta = 0.1/1000;
d = 1; S = pi()*d^2/4;
Q = 4800/3600;
L = 100*10^3;
%parametry nasosa
a0 = 900; b = 5*10^(-6); n_nasosa = 50;
w0 = 3000; w = 3000;

c = 1000; %skorost' zvuka
N = 100; dx = L/N; dt = dx/c;
x = 0:dx:L;

v0 = zeros(1,N+1); v0(:) = 4*Q/pi()/d^2;
p0 = zeros(1,N+1);

for i = 1:N+1
   p0(i) =  p_nach - (p_nach-p_kon)/ L * dx*(i-1);
end
p_vs = p0(n_nasosa); p_nag = p0(n_nasosa);

for j = 1:1000
t = t + dt;
if t > 100
    w = 3000 - 2000/10*(t-100);
end
if t > 110
    w = 1000;
end
if t > 400
    w = 1000 + 2000/10*(t-400);
end
if t > 410
    w = 3000;
end
a = a0*(w/w0)^2;
% if t > 100
%     p_kon = 5*10^6;
% end
% if t > 120
%     p_nach = 3*10^6;
% end
p1 = zeros(1,N+1); v1 = zeros(1,N+1);
%levaya granica
i = 1;
p1(i) = p_nach;
InvB = p0(i+1) - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
v1(i) = (p1(i) - InvB) / (rho*c);

%tsentralnaya chast truby do nasosa
for i = 2:n_nasosa-2
   InvA = p0(i-1) + rho*c*v0(i-1) - lambda(v0(i-1),d,visc,delta)*rho*v0(i-1)*abs(v0(i-1))/2/d*c*dt;
   InvB = p0(i+1) - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
   p1(i) = (InvA+InvB)/2;
   v1(i) = (InvA-InvB)/(2*rho*c);
end
%pered nasosom
   i = n_nasosa-1;
   InvA = p0(i-1) + rho*c*v0(i-1) - lambda(v0(i-1),d,visc,delta)*rho*v0(i-1)*abs(v0(i-1))/2/d*c*dt;
   InvB = p_vs - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
   p1(i) = (InvA+InvB)/2;
   v1(i) = (InvA-InvB)/(2*rho*c);
%nasos
i = n_nasosa;
InvA = p0(i-1) + rho*c*v0(i-1) - lambda(v0(i-1),d,visc,delta)*rho*v0(i-1)*abs(v0(i-1))/2/d*c*dt;
InvB = p0(i+1) - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
discr = 4*(c/9.81)^2 - 4*((InvB-InvA)/rho/9.81-a)*b*S^2*3600^2;
v1(i) = (-2*c/9.81 + sqrt(discr)) / (2*b*S^2*3600^2);
p_vs = InvA - rho*c*v1(i);
p_nag = InvB +rho*c*v1(i);
p1(n_nasosa) = p_vs;
%secheniye posle nasosa
   i = n_nasosa+1;
   InvA = p_nag + rho*c*v0(i-1) - lambda(v0(i-1),d,visc,delta)*rho*v0(i-1)*abs(v0(i-1))/2/d*c*dt;
   InvB = p0(i+1) - rho*c*v0(i+1) + lambda(v0(i+1),d,visc,delta)*rho*v0(i+1)*abs(v0(i+1))/2/d*c*dt;
   p1(i) = (InvA+InvB)/2;
   v1(i) = (InvA-InvB)/(2*rho*c);

   %truba posle nasosa
for i = n_nasosa+2:N
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