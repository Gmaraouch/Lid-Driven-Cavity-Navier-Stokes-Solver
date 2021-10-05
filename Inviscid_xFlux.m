function F=Inviscid_xFlux(rho, u, v, E,p)

F=cell(4,1);
F{1}=rho.*u;
F{2}=rho.*u.*u+p;
F{3}=rho.*u.*v;
F{4}=u.*(rho.*E+p);

end