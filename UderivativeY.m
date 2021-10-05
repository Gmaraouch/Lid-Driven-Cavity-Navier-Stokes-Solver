function [dudy]=UderivativeY(u,dy,Periodic)
    
    dudy=u;

    if(strcmp('yes',Periodic))
        
        for i=1:size(u,2)

            if(i==size(u,2))
            dudy(i,:)=DudxCentral2(u(1,:),u(end-1,:),dy);
            elseif(i==1)
            dudy(i,:)=DudxCentral2(u(2,:),u(end,:),dy);
            else
            dudy(i,:)=DudxCentral2(u(i+1,:),u(i-1,:),dy);
            end

        end
        
    else
        
        for i=1:size(u,2)

            if(i==size(u,2))
            %dudy(i,:)=(u(end-2,:)-4*u(end-1,:)+3*u(end,:))/dy;
            dudy(i,:)=DudxBackward2(u(end-2,:),u(end-1,:),u(end,:),dy);
            elseif(i==1)
            %dudy(i,:)=(-u(3,:)+4*u(2,:)-3*u(1,:))/dy;
            dudy(i,:)=DudxForward2(u(3,:),u(2,:),u(1,:),dy);
            else
            %dudy(i,:)=(u(i+1,:)-u(i-1,:))/(2*dy);
            dudy(i,:)=DudxCentral2(u(i+1,:),u(i-1,:),dy);
            end

        end
    end
    
end