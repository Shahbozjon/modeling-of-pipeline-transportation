clc; clear
%raschet slozhnoy truboprovodnoy seti

L = [15,24,20,32,65,32,52]*10^3;
d = [1,1,0.8,1,0.7,1,0.6];
z = [100,100,100,100,100,100,100,100];
delta = 0.1/1000;
rho = 800; visc = 12/10^6; 

dp = 1*10^3; %popravka
p = [20,0,0,0,65,19,12,36] * 10^5;
Q = zeros(1,7);
%-----------konets ishodnih dannyh

%step1
p(2) = 50*10^5; p(3) = 25*10^5; p(4) = 35*10^5;

dQ2 = 1; dQ3 = 1; dQ4 = 1; %dummy line
i = 0; %iterator
while (abs(dQ2) + abs(dQ3) + abs(dQ4)) > 0.01
i = i+1;
%step2
Q(1) = find_Q(L(1),d(1),rho,visc,delta,p(1),p(2),z(1),z(2));
Q(2) = find_Q(L(2),d(2),rho,visc,delta,p(2),p(3),z(2),z(3));
Q(3) = find_Q(L(3),d(3),rho,visc,delta,p(2),p(4),z(2),z(4));
Q(4) = find_Q(L(4),d(4),rho,visc,delta,p(3),p(5),z(3),z(5));
Q(5) = find_Q(L(5),d(5),rho,visc,delta,p(3),p(6),z(3),z(6));
Q(6) = find_Q(L(6),d(6),rho,visc,delta,p(4),p(7),z(4),z(7));
Q(7) = find_Q(L(7),d(7),rho,visc,delta,p(4),p(8),z(4),z(8));

%step3
dQ2 = Q(1)-Q(2)-Q(3);
dQ3 = Q(2)-Q(4)-Q(5);
dQ4 = Q(3)-Q(6)-Q(7);

%step5
if dQ2 > 0
   p(2) = p(2) + dp; 
else
   p(2) = p(2) - dp;
end

if dQ3 > 0
   p(3) = p(3) + dp; 
else
   p(3) = p(3) - dp;
end

if dQ4 > 0
   p(4) = p(4) + dp; 
else
   p(4) = p(4) - dp;
end

end

disp([dQ2 dQ3 dQ4])
disp(i)
