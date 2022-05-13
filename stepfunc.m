%step function which decreases from max to min at a certain onset value

function stepfunction=stepfunc(tonset,timevec,minv,maxv)
    tspan=length(timevec);
    for i=1:tspan
        if (timevec(i))>tonset
            stepfunction(i)=minv;
        else
            stepfunction(i)=maxv;
        end
    end 
end