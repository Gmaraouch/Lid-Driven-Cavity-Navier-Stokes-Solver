function [p]=Update_Pressure(w,Constant)

gamma=Constant.gamma;

rho=w{1};
rhou=w{2};
rhov=w{3};
rhoE=w{4};

u=rhou./rho;
v=rhov./rho;

p=(gamma-1)*(rhoE-(1/2).*rho.*(u.^2+v.^2));

end