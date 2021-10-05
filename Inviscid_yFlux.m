function G=Inviscid_yFlux(rho, u, v, E,p)

G=cell(4,1);

G{1}=rho.*v;
G{2}=rho.*v.*u;
G{3}=rho.*v.*v+p;
G{4}=v.*(rho.*E+p);

end