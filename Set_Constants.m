function Constant=Set_Constants(dx,dy, Nx, Ny)

gamma=1.4;
Ma=0.1;
Lx=1;
Ly=1;
Re=100;
Pr=0.71;
Rs=1;
Uw=1;
dt=0.001;
Initial_density=1;
To=(Uw/Ma)^2/gamma/Rs;
nu=abs(Uw)*Lx/Re;
alpha=nu/Pr;
mu=nu*Initial_density;
cp=(gamma)/(gamma-1)*Rs;
k=cp*mu/Pr;

Constant.gamma=gamma;
Constant.Ma=Ma;
Constant.Lx=Lx;
Constant.Ly=Ly;
Constant.Re=Re;
Constant.Pr=Pr;
Constant.Rs=Rs;
Constant.Uw=Uw;
Constant.mu=mu;
Constant.dt=dt;
Constant.dx=dx;
Constant.dy=dy;
Constant.Nx=Nx;
Constant.Ny=Ny;
Constant.To=To;
Constant.nu=nu;
Constant.alpha=alpha;
Constant.k=k;
Constant.cp=cp;
Constant.Initial_density=Initial_density;

end