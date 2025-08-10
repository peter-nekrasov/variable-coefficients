function [r1,grad,hess,third] = planewave(k,r,theta)

k1 = k*cos(theta);
k2 = k*sin(theta);

r1 = exp(1i*(k1*r(1,:)+k2*r(2,:))).';

gx = r1.*(1i*k1);
gy = r1.*(1i*k2);
grad = [gx, gy];

hxx = r1.*(1i*k1).^2; 
hxy = r1.*(1i*k1)*(1i*k2); 
hyy = r1.*(1i*k2).^2;
hess = [hxx, hxy, hyy];

txxx = r1.*(1i*k1).^3;
txxy = r1.*(1i*k1).^2*(1i*k2); 
txyy = r1.*(1i*k1)*(1i*k2).^2; 
tyyy = r1.*(1i*k2).^3; 
third = [txxx, txxy, txyy, tyyy];

end