function Q = find_Q(L,d,rho,visc,delta,p0,pL,z0,zL)
%script opredelyaet rashod nefti v truboprovode
%s zadannimy parametrami
g = 9.81;

lambda_old = 0.02;
lambda_new = 0.01; i =1;
while (abs(lambda_old-lambda_new)/lambda_old > 0.01)
i = i+1;
v = sign(p0-pL) * sqrt(((abs(p0-pL))/rho/g + z0 - zL) * 2*d*g/lambda_old/L);
lambda_new = lambda(v,d,visc,delta);
lambda_old = lambda_new;
end
v_result = v;
Q = v_result * pi() * d^2/4;

end