function [wt,Finv,Ginv,Fvis,Gvis]=fRK44(w,Constant,Periodic,Enforce_Boundaries,Viscous)
    

    dx=Constant.dx;
    dy=Constant.dy;
    dt=Constant.dt;
    
    w2=cell(4,1);
    w3=w2;
    w4=w2;
    wt=w2;
    Rw1=cell(4,1);
    Rw2=Rw1;
    Rw3=Rw1;
    Rw4=Rw1;
    
    w1=w;
    if(strcmp(Enforce_Boundaries,'yes'))
    [w1]=Set_Boundaries(w1,Constant);
    end
    
    if(strcmp('yes',Viscous))
        [F1inv,F1vis]=xflux(w1,Constant,Viscous,Periodic,Enforce_Boundaries);
        [G1inv,G1vis]=yflux(w1,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw1inv=R_Central_2(F1inv,G1inv,dx,dy,Periodic);
        Rw1vis=R_Central_2(F1vis,G1vis,dx,dy,Periodic);

        for i=1:4
        Rw1{i}=Rw1inv{i}+Rw1vis{i};  
        end
    else
        F1inv=xflux(w1,Constant,Viscous,Periodic,Enforce_Boundaries);
        G1inv=yflux(w1,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw1inv=R_Central_2(F1inv,G1inv,dx,dy,Periodic);

        Rw1=Rw1inv;  
    end
    
    for i=1:4
    w2{i}=w{i}+(dt/2).*Rw1{i};
    end
    
    if(strcmp(Enforce_Boundaries,'yes'))
    [w2]=Set_Boundaries(w2,Constant);
    end
    
    if(strcmp('yes',Viscous))
        
        [F2inv,F2vis]=xflux(w2,Constant,Viscous,Periodic,Enforce_Boundaries);
        [G2inv,G2vis]=yflux(w2,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw2inv=R_Central_2(F2inv,G2inv,dx,dy,Periodic);
        Rw2vis=R_Central_2(F2vis,G2vis,dx,dy,Periodic);

        for i=1:4
        Rw2{i}=Rw2inv{i}+Rw2vis{i}; 
        end
    
    else
        
        F2inv=xflux(w2,Constant,Viscous,Periodic,Enforce_Boundaries);
        G2inv=yflux(w2,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw2inv=R_Central_2(F2inv,G2inv,dx,dy,Periodic);

        Rw2=Rw2inv;
        
    end
        
    for i=1:4
    w3{i}=w{i}+(dt/2).*Rw2{i};
    end
    
    if(strcmp(Enforce_Boundaries,'yes'))
    [w3]=Set_Boundaries(w3,Constant);
    end
    
    if(strcmp('yes',Viscous))
        
        [F3inv,F3vis]=xflux(w3,Constant,Viscous,Periodic,Enforce_Boundaries);
        [G3inv,G3vis]=yflux(w3,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw3inv=R_Central_2(F3inv,G3inv,dx,dy,Periodic);
        Rw3vis=R_Central_2(F3vis,G3vis,dx,dy,Periodic);

        for i=1:4
        Rw3{i}=Rw3inv{i}+Rw3vis{i};
        end
    
    else

    F3inv=xflux(w3,Constant,Viscous,Periodic,Enforce_Boundaries);
    G3inv=yflux(w3,Constant,Viscous,Periodic,Enforce_Boundaries);
    Rw3inv=R_Central_2(F3inv,G3inv,dx,dy,Periodic);
    
    Rw3=Rw3inv;
    
    end
    
    for i=1:4
    w4{i}=w{i}+(dt).*Rw3{i};
    end
    
    if(strcmp(Enforce_Boundaries,'yes'))
    [w4]=Set_Boundaries(w4,Constant);
    end
    
    if(strcmp('yes',Viscous))
        [F4inv,F4vis]=xflux(w4,Constant,Viscous,Periodic,Enforce_Boundaries);
        [G4inv,G4vis]=yflux(w4,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw4inv=R_Central_2(F4inv,G4inv,dx,dy,Periodic);
        Rw4vis=R_Central_2(F4vis,G4vis,dx,dy,Periodic);

        for i=1:4
        Rw4{i}=Rw4inv{i}+Rw4vis{i};
        end
    else
        [F4inv]=xflux(w4,Constant,Viscous,Periodic,Enforce_Boundaries);
        [G4inv]=yflux(w4,Constant,Viscous,Periodic,Enforce_Boundaries);
        Rw4inv=R_Central_2(F4inv,G4inv,dx,dy,Periodic);

        Rw4=Rw4inv;
    end
    
    for i=1:4
    wt{i}=w{i}+(dt/6).*(Rw1{i}+2.*Rw2{i}+2.*Rw3{i}+Rw4{i});
    end
    
    if(strcmp(Enforce_Boundaries,'yes'))
    [wt]=Set_Boundaries(wt,Constant);
    end
    
    if(strcmp('yes',Viscous))
    [Finv,Fvis]=xflux(w,Constant,Viscous,Periodic,Enforce_Boundaries);
    [Ginv,Gvis]=yflux(w,Constant,Viscous,Periodic,Enforce_Boundaries);
    else
    [Finv]=xflux(w,Constant,Viscous,Periodic,Enforce_Boundaries);
    [Ginv]=yflux(w,Constant,Viscous,Periodic,Enforce_Boundaries);        
    end
end