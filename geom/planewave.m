function pw = planewave(k,r,theta,nterm)

if nargin < 4
    nterm = 1;
end

k1 = k*cos(theta);
k2 = k*sin(theta);

r1 = exp(1i*(k1*r(1,:)+k2*r(2,:))).';

gx = r1.*(1i*k1);
gy = r1.*(1i*k2);

hxx = r1.*(1i*k1).^2; 
hxy = r1.*(1i*k1)*(1i*k2); 
hyy = r1.*(1i*k2).^2;

txxx = r1.*(1i*k1).^3;
txxy = r1.*(1i*k1).^2*(1i*k2); 
txyy = r1.*(1i*k1)*(1i*k2).^2; 
tyyy = r1.*(1i*k2).^3; 

pw = cat(3,r1,gx,gy,hxx,hxy,hyy,txxx,txxy,txyy,tyyy);
pw = pw(:,:,1:nterm);

end