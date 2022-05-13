
%function defining how the K fraction changes with time after onset of K
%change 


function kfraction=kfunc(tkonset,timep,Kmin,Kmax,thalf,h)
    tspan1=length(timep);
    if tspan1==1
        kfraction=Kmin+(Kmax*(timep-tkonset)^h)/((thalf-tkonset)^h+(timep-tkonset)^h);
    else
        for j=1:tspan1
            if (timep(j))>tkonset
                kfraction(j)=(Kmin+(Kmax*(timep(j)-tkonset)^h)/((thalf-tkonset)^h+(timep(j)-tkonset)^h));
            else
                kfraction(j)=Kmin;
            end
        end 
    end 
end 