function [u,v,V]=Get_Velocity(w)

rho=w{1};
rhou=w{2};
rhov=w{3};

u=rhou./rho;
v=rhov./rho;
V=sqrt(u.^2+v.^2);


end