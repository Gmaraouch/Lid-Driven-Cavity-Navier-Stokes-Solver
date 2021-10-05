function [w]=Isentropic_Vortex_Initial(Constant,x,y)
    
    xc=Constant.xc;
    yc=Constant.yc;
    S=Constant.S;
    gamma=Constant.gamma;
    Ma=Constant.Ma;
    R=Constant.R;
    
    f=(1-(x-xc).^2-(y-yc).^2)/(2*R^2);
    rho=(1-((S^2*Ma^2*(gamma-1)).*exp(2*f))/(8*pi^2)).^(1/(gamma-1));
    u=S*(y-yc).*exp(f)/(2*pi*R);
    v=1-S*(x-xc).*exp(f)/(2*pi*R);
    p=(rho.^(gamma))./(gamma*Ma^2);
    
    rhoE=p./(gamma-1)+1/2*rho.*(u.^2+v.^2);
    
    w=cell(4,1);
    w{1}=rho;
    w{2}=rho.*u;
    w{3}=rho.*v;
    w{4}=rhoE;
    
end