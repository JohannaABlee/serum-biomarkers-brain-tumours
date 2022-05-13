%function to calculate the volume of necrositic region assuming necrotic
%region grows via gompertzian 
function volumen=voln(tnonset,timepoint,maxvoln,an,Vn0)
    tspan=length(timepoint);
    for i=1:tspan
        if (timepoint(i))>tnonset
            volumen=maxvoln*exp(log(Vn0/maxvoln)*exp(-an*(timepoint-tnonset)));%necrosis volume can be expressed as number of cells or volume
        else
            volumen=0;
        end
    end 
end

