function prodn=Un(qn,rn,tnonset,timepoint,maxvoln,an,Vn0)
    tspan=length(timepoint);
    for i=1:tspan
        if (timepoint(i))>tnonset
            volumen=maxvoln*exp(log(Vn0/maxvoln)*exp(-an*(timepoint-tnonset)));%necrosis volume can be expressed as number of cells or volume
            prodn=qn*rn*volumen*log(maxvoln/volumen);
        else
            volumen=0;
            prodn=0;
        end
       
    end 
end