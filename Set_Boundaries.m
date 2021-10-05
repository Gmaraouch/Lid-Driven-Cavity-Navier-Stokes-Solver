function [w,p,T]=Set_Boundaries(w,Constant)
%set rhou and rhov at the wall to zero

w{2}(:,1)=0;
w{2}(1,:)=0;
w{2}(:,end)=0;
if(strcmp('no',Constant.No_Slip_At_Corners))
w{2}(end,1:end)=Constant.Uw*w{1}(end,1:end);
else
w{2}(end,2:end-1)=Constant.Uw*w{1}(end,2:end-1);
end

w{3}(:,1)=0;
w{3}(1,:)=0;
w{3}(:,end)=0;
w{3}(end,:)=0;

p=Update_Pressure(w,Constant);

p(:,1)=p(:,2);
p(1,:)=p(2,:);
p(:,end)=p(:,end-1);
p(end,:)=p(end-1,:);

w{1}(:,1)=p(:,2)/(Constant.Rs*Constant.To);
w{1}(1,:)=p(2,:)/(Constant.Rs*Constant.To);
w{1}(:,end)=p(:,end-1)/(Constant.Rs*Constant.To);
w{1}(end,:)=p(end-1,:)/(Constant.Rs*Constant.To);

T=p./w{1}./Constant.Rs;
T(:,1)=Constant.To;
T(1,:)=Constant.To;
T(:,end)=Constant.To;
T(end,:)=Constant.To;

rho=w{1};
rhou=w{2};
rhov=w{3};

u=rhou./rho;
v=rhov./rho;

rhoE=p./(Constant.gamma-1)+1/2*rho.*(u.^2+v.^2);

w{4}=rhoE;

end