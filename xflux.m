%This function calculates the flux in the x direction.
%It returns the Inviscid Flux and the Viscous Flux

function [Finv,Fv]=xflux(w,Constant,Viscous,Periodic,Enforce_Boundaries)

%w=[rho, rho*u,rho*v,rho*E];
%F=[rho*u,rho*u^2+p,rho*u*v,u*(rho*E+p)];

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

Finv=Inviscid_xFlux(rho, u, v, E, p);

if(strcmp('yes',Viscous))
Fv=Viscous_xFlux(u, v, T, Constant,Periodic);
end

end