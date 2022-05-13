% function for tumour growth assuming grompertzian mode of growth
%input parameters are the time point, the max volume(maxvol), growth
%rate(a) and strating volume V0.

function volumeT=volT(timepoint,maxvol,a,V0)
    volumeT=maxvol*exp(log(V0/maxvol)*exp(-a*timepoint));%tumour volume
end