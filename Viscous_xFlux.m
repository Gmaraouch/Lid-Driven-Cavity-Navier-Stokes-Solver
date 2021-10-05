function F=Viscous_xFlux(u, v, T,Constant,Periodic)

mu=Constant.mu;
dx=Constant.dx;
dy=Constant.dy;

dudx=UderivativeX(u,dx,Periodic);
dudy=UderivativeY(u,dy,Periodic);

dvdx=UderivativeX(v,dx,Periodic);
dvdy=UderivativeY(v,dy,Periodic);

Tauxx=2/3*mu.*(2*dudx-dvdy);
Tauyy=2/3*mu.*(2*dvdy-dudx);
Tauxy=mu.*(dudy+dvdx);
Tauyx=Tauxy;

qx=-Constant.k*UderivativeX(T,dx,Periodic);

F=cell(4,1);
F{1}=zeros(size(u,1),size(u,2));
F{2}=-Tauxx;
F{3}=-Tauxy;
F{4}=qx-u.*Tauxx-v.*Tauyx;

end