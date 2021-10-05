function [w]=Initial_Variables(Constant,x,y)

    w=cell(4,1);
    Nx=length(x);
    Ny=length(y);
    rho(1:Nx,1:Ny)=Constant.Initial_density;
    u(1:Nx,1:Ny)=0;
    
    if(strcmp('no',Constant.No_Slip_At_Corners))
    u(end,1:end)=Constant.Uw;
    else
    u(end,2:end-1)=Constant.Uw;
    end
    v(1:Nx,1:Ny)=0;
    
    p=rho.*Constant.Rs.*Constant.To;
   
    rhoE=p./(Constant.gamma-1)+1/2*rho.*(u.^2+v.^2);

    w{1}=rho;
    w{2}=rho.*u;
    w{3}=rho.*v;
    w{4}=rhoE;

    %[w]=Set_Boundaries(w,Constant);
end