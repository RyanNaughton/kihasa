
function CES=CES(h,minv,tinv,gamma,rho,phi,eps_h)

CES=(gamma(1)*h.^phi + gamma(2)*minv.^phi + (1-gamma(1)-gamma(2))*tinv.^phi).^(rho/phi)+eps_h;

%CES=(gamma(1)*h.^phi + (1-gamma(1))*minv.^phi + 0*tinv.^phi).^(rho/phi)+eps_h;