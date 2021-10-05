function G=Viscous_yFlux(u, v, T,Constant,Periodic)

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

qy=-Constant.k*UderivativeY(T,dy,Periodic);

G=cell(4,1);
G{1}=zeros(size(u,1),size(u,2));
G{2}=-Tauxy;
G{3}=-Tauyy;
G{4}=qy-u.*Tauxy-v.*Tauyy;

end