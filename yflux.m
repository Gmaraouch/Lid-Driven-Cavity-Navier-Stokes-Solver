%This function calculates the flux in the y direction.
%It returns the Inviscid Flux and the Viscous Flux

function [Ginv,Gv]=yflux(w,Constant,Viscous,Periodic,Enforce_Boundaries)

%w=[rho, rho*u,rho*v,rho*E];
%G=[rho*v,rho*u*v,rho*v^2+p,v*(rho*E+p)];
if(strcmp('yes',Enforce_Boundaries))
[w,p,T]=Set_Boundaries(w,Constant);
else
p=Update_Pressure(w,Constant);    
end

rho=w{1};
rhou=w{2};
rhov=w{3};
rhoE=w{4};

u=rhou./rho;
v=rhov./rho;
E=rhoE./rho;


Ginv=Inviscid_yFlux(rho, u, v, E,p);

if(strcmp('yes',Viscous))
Gv=Viscous_yFlux(u, v, T, Constant,Periodic);
end

end