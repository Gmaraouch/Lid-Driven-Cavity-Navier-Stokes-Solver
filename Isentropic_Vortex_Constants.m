function [Constant]=Isentropic_Vortex_Constants(dx,dy,Nx,Ny)

Lx=20;
Ly=20;
gamma=1.4;
R=1.5;
Ma=0.4;
xc=10;
yc=10;
S=13.5;
dt=0.01;

Constant.gamma=gamma;
Constant.Ma=Ma;
Constant.Lx=Lx;
Constant.Ly=Ly;
Constant.dt=dt;
Constant.dx=dx;
Constant.dy=dy;
Constant.Nx=Nx;
Constant.Ny=Ny;
Constant.S=S;
Constant.xc=xc;
Constant.yc=yc;
Constant.R=R;
    
end