%%This code runs a flow in a cavity with a moving plate on the top

clear 
clc
clf

%addpath('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 4\Question 1');

%Note that you should write 'yes' with no capitals
%Define if the problem has periodic boundaries or no
Periodic='no';
%Define if the problem has viscosity in it or not (i.e. Euler or
%Navier-Stokes)
Viscous='yes';
%Do you have any boundary conditions that need to be held?
Enforce_Boundaries='yes';
%Save each iterated solution? Be careful, if the grid spacing is too big,
%the computer will run out of RAM
Save_Time='no';

%Define the boundaries length
Lx=1;
Ly=1;

%Define the number of grid spacing in the x-direction (assume domain is
%square)
Nx=[33,65,129];
%%%%%%
%If you would like to run the isentropic vortex example, change the string
%to 'yes' (cap sensitive)
Isentropic_vortex='no';
if(strcmp('yes',Isentropic_vortex))
    Periodic='yes';
    Viscous='no';
    Enforce_Boundaries='no';
    clear Nx
    Nx=[60,120,240,480];
    Lx=20;   
    Ly=20;
end
%%%%

Ny=Nx;

for i=1
    
    %Calculate the grid spacing and define the domain x and y
    if(strcmp('yes',Periodic))
    dx=Lx/Nx(i);
    dy=Ly/Ny(i);
    x=0:dx:(Lx-dx);
    y=0:dy:(Ly-dx);
    else
    dx=Lx/(Nx(i)-1);
    dy=Ly/(Ny(i)-1);
    x=0:dx:(Lx);
    y=0:dy:(Ly);
    end
    
    %Set the constants for the problem
    if(strcmp(Isentropic_vortex,'no'))
    Constant=Set_Constants(dx,dy,Nx(i),Ny(i));
    else
    Constant=Isentropic_Vortex_Constants(dx,dy,Nx(i),Ny(i));
    end
    
    %Is the no-slip condition held at the corners where the moving plate
    %is? If yes, change the string to 'yes'.
    Constant.No_Slip_At_Corners='yes';
    
    %Time of the simulation
    time=20;
    
    %Calculate the number of increments of the simulation
    nt=cast(time/Constant.dt,'int16');
    
    %Create a mesh of the grid
    [x,y]=meshgrid(x,y);
    
    %Calculate the initial state properties (rho, rhou, rhov and rhoE)
    if(strcmp('no',Isentropic_vortex))
    wo=Initial_Variables(Constant,x,y);
    else
    wo=Isentropic_Vortex_Initial(Constant,x,y);
    end

    
    %Calculate the initial flux in the x and y direction.
    if(strcmp('no',Viscous))
        
    [Finv]=xflux(wo,Constant,Viscous,Periodic,Enforce_Boundaries);
    [Ginv]=yflux(wo,Constant,Viscous,Periodic,Enforce_Boundaries);
    
    wt=cell(nt,1);
    wt{1}.w=wo;
    wt{1}.Finv=Finv;
    wt{1}.Ginv=Ginv;
    
    else
    [Finv,Fvis]=xflux(wo,Constant,Viscous,Periodic,Enforce_Boundaries);
    [Ginv,Gvis]=yflux(wo,Constant,Viscous,Periodic,Enforce_Boundaries);
    
    wt=cell(nt,1);
    wt{1}.w=wo;
    wt{1}.Finv=Finv;
    wt{1}.Ginv=Ginv;
    wt{1}.Fvis=Fvis;
    wt{1}.Gvis=Gvis; 
    end
    
    w=wo;
    %Run the simulation
    if(strcmp('yes',Viscous))
        for t=1:nt
            if(strcmp('yes',Save_Time))
                [wt{t+1}.w,wt{t+1}.Finv,wt{t+1}.Ginv,wt{t+1}.Fvis,wt{t+1}.Gvis]=fRK44(w,Constant,Periodic,Enforce_Boundaries,Viscous);
                [w,p,T]=Set_Boundaries(wt{t+1}.w,Constant);
                drawnow;
                contourf(x,y,wt{t+1}.w{1});
                %contourf(x,y,wt{t+1}.w{2}./wt{t+1}.w{1});
                %surf(x,y,wt{t+1}.w{1},'EdgeColor','None');
                %w=wt{t+1}.w;
                %pause(2)
                
                %save the density as a csv file
                %s=strcat('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 4\Question 1\2nd Order\Density\',...
                       %  'Nx=', num2str(Nx(i)),'\Density ',num2str(t),'.csv');
                %csvwrite(s,wt{t+1}.w{1});

            else
                [wt,Finv,Ginv,Fvis,Gvis]=fRK44(w,Constant,Periodic,Enforce_Boundaries,Viscous);
                [w,p,T]=Set_Boundaries(wt,Constant);
                drawnow;
                contourf(x,y,w{1});
                %surf(x,y,wt{1},'EdgeColor','None');
                %w=wt;
            end
        end
    else
        for t=1:nt
            if(strcmp('yes',Save_Time))
                [wt{t+1}.w,wt{t+1}.Finv,wt{t+1}.Ginv]=fRK44(w,Constant,Periodic,Enforce_Boundaries,Viscous);
                drawnow;
                %contourf(x,y,wt{t+1}.w{1});
                surf(x,y,wt{t+1}.w{1},'EdgeColor','None');
                w=wt{t+1}.w;
            else
                %[wt,Finv,Ginv]=fRK44_Inviscid(w,Constant,Periodic,Enforce_Boundaries);
                [wt,Finv,Ginv]=fRK44(w,Constant,Periodic,Enforce_Boundaries,Viscous);
                drawnow;
                surf(x,y,wt{1},'EdgeColor','None');
                w=wt;
            end
                
        end
    end
    
    figure(1)
    contourf(x,y,p);
    xlabel('x (units)');
    ylabel('y (units)');
    title(['Pressure Contour for a Grid Spacing of ' num2str(Nx(i)) 'x' num2str(Ny(i))]);
    colorbar;
    print('-f1',strcat('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 4\Question 1\2nd Order\Plots\'...
            ,'Pressure Contour for Nx=',num2str(Nx(i))),'-dtiff','-painters');
    
    [u,v,V]=Get_Velocity(w);
    contourf(x,y,V);
    xlabel('x (units)');
    ylabel('y (units)');
    title(['Velocity Magnitude Contour for a Grid Spacing of ' num2str(Nx(i)) 'x' num2str(Ny(i))]);
    colorbar;  
    print('-f1',strcat('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 4\Question 1\2nd Order\Plots\'...
            ,'Velocity Magnitude Contour for Nx=',num2str(Nx(i))),'-dtiff','-painters');
    
    %values for the vertical velocity over the x-range taken from Ghia et.
    %al (1982)
    xt=[1.0000;0.9688;0.9609;0.9531;0.9453;0.9063;0.8594;0.8047;0.5000;0.2344;0.2266;0.1563;0.0938;0.0781;0.0703;0.0625;0.0000];    
    vt=[0.00000;-0.05906;-0.07391;-0.08864;-0.10313;-0.16914;-0.22445;-0.24533;0.05454;0.17527;0.17507;0.16077;0.12317;0.10890;0.10091;0.09233;0.00000];    
    
    half=cast(Ny(i)/2,'int16');
    plot(xt,vt,x(1,:),v(half,:));
    xlabel('x (units)');
    ylabel('Vertical Velocity Component (units)');
    title(['Vertical Velocity of the Lid Driving Cavity at the Center with Nx=' num2str(Nx(i))]);
    legend('Reference Vertical Velocity',['Simulated Vertical Velocity with Nx=' num2str(Nx(i))]);
    print('-f1',strcat('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 4\Question 1\2nd Order\Plots\'...
            ,'Vertical Velocity for Nx=',num2str(Nx(i))),'-dtiff','-painters');
        
    yt=[1.0000;0.9766;0.9688;0.9609;0.9531;0.8516;0.7344;0.6172;0.5000;0.4531;0.2813;0.1719;0.1016;0.0703;0.0625;0.0547;0.0000];
    ut=[1.00000;0.84123;0.78871;0.73722;0.68717;0.23151;0.00332;-0.13641;-0.20581;-0.21090;-0.15662;-0.10150;-0.06434;-0.04775;-0.04192;-0.03717;0.00000];

    plot(yt,ut,y(:,1),u(:,half));
    xlabel('y (units)');
    ylabel('Horizontal Velocity Component (units)');
    title(['Horizontal Velocity of the Lid Driving Cavity at the Center with Nx=' num2str(Nx(i))]);
    legend('Reference Hozitonal Velocity',['Simulated Horizontal Velocity with Nx=' num2str(Nx(i))],'Location','Northwest');
    print('-f1',strcat('C:\Users\Gmara\OneDrive - Concordia University - Canada\Concordia\Fall 2019\ENGR 6251\Assignment 4\Question 1\2nd Order\Plots\'...
            ,'Horizontal Velocity for Nx=',num2str(Nx(i))),'-dtiff','-painters');
end

%{
V=cell(length(wt),1);
for i=1:length(wt)
    
    [V{i}.u,V{i}.v]=Get_Velocity(wt{i}.w);

end

for i=1:length(wt)

drawnow;
quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),V{i}.u(1:2:end,1:2:end),V{i}.v(1:2:end,1:2:end),2);
xlim([0 1]);
ylim([0 1]);
end
%}