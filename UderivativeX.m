function [dudx]=UderivativeX(u,dx,Periodic)
    
    dudx=u;

    if(strcmp('yes',Periodic))
        
        for i=1:size(u,2)

            if(i==1)
            dudx(:,i)=DudxCentral2(u(:,2),u(:,end),dx);
            elseif(i==size(u,2))
            dudx(:,i)=DudxCentral2(u(:,1),u(:,end-1),dx);
            else
            dudx(:,i)=DudxCentral2(u(:,i+1),u(:,i-1),dx);
            end

        end
        
    else
        for i=1:size(u,2)

            if(i==1)
            %dudx(:,i)=(-u(:,3)+4*u(:,2)-3*u(:,1))/dx;
            dudx(:,i)=DudxForward2(u(:,3),u(:,2),u(:,1),dx);
            elseif(i==size(u,2))
            %dudx(:,i)=(u(:,end-2)-4*u(:,end-1)+3*u(:,end))/dx;
            dudx(:,i)=DudxBackward2(u(:,end-2),u(:,end-1),u(:,end),dx);
            else
            %dudx(:,i)=(u(:,i+1)-u(:,i-1))/(2*dx);
            dudx(:,i)=DudxCentral2(u(:,i+1),u(:,i-1),dx);
            end

        end
    end
end