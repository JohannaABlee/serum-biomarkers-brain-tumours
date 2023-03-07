%Code for mathematical modelling of serum GFAP with glioma growth 
%Johanna Blee 04/04/2021

%-----parameter intialisation ----
time=linspace(1,600,600);%time of growth simulation
pt=5714*10^4;%tumour cell density (cells/ml)
pn=4800*10^4;%necrotic cell density (cells/ml)
N0=1; %intial number of tumour cells 
V0=1/pt; %intial volume of tumour cells
Nn0=1;
Vn0=1/pn;
Vp=4500; %ml mean plasma volume
cthr=0.12;%current ut-off limit
CPss=0.012; %ng/ml average conc in serum of healthy patients
yp=0.7;%days-1 decay rate in serum
KHUH=CPss*(yp*Vp); %average production*fraction in healthy patients
qt0=CPss*Vp%initial value in healthy induvidual at tumour onset

%growth of tumour and necrosis
a=0.008;%growth rate of tumour 
b=0.009;%growth rate of necrotic regions
maxVt=158;%maximum tumour volume (ml)
maxNt=maxVt*pt; %maximum number of tumour cells
maxVn=150;%maximum necrotic volume
maxNn=pn*maxVn;%maximum number of necrotic cells

Vt=volT(time,maxVt,a,V0); %tumour vol ml 
Nt=volT(time,maxNt,a,N0); %number of tumour cells

vn0=2;% volume of onset of necoris ml (0.5cm^3);
nn0=vn0*pt; %number of tumour cells at onset of necrosis 
[minValue,closestIndex]=min(abs(Vt-vn0));
tn0=time(closestIndex); %time onset necrosis

Vn=voln(tn0,time,maxVn,b,Vn0); %necrotic volume
Nn=voln(tn0,time,maxNn,b,Nn0); %number of necrotic cells


%necrosis constants
Rd=b;%rate of cell death- same as necrosis growth rate
Qn=3.1*10^-4;%average level of GFAP in necrotic cells

%kfunction intilaisation
tk0=tn0%inthis scenario assume necrosis and k onset align
Ktmin=0;%minimum k value
Ktmax=0.5;%maximum k value
h=7;%rate of k change hill slope constant
Ktt=225;%half time - when k is half max
kfunction=kfunc(tk0,time,Ktmin,Ktmax,Ktt,h);

%----ODE to give GFAP serum levels as a function of time ----
[time2,ytn1]=ode45(@(time2,ytn1)((kfunc(tk0,time2,Ktmin,Ktmax,Ktt,h)*Un(Qn,Rd,tn0,time2,maxNn,b,Nn0))+KHUH-ytn1*yp),time,qt0);
CP1=ytn1/Vp; %concentration of GFAP with time
vol2=volT(time2,maxVt,a,V0); %tumour volume at a given time point


%---time and voloume when critical threshold reached----
[minValue2,closestIndex2] = min(abs(CP1-cthr)); 
tthresh1=time2(closestIndex2); %time threshold reached
vthresh1=volT(tthresh1,maxVt,a,V0); %detection volume tumour volume at which current GFAP threshold is reached 


%-- Figures for GFAP with volume and with time ----
figure
plot(vol2,CP1)
ylabel('Serum GFAP (ng/ml)')
xlabel('Tumour volume (ml)')


figure
plot(time2,CP1)
ylabel('Serum GFAP (ng/ml)')
xlabel('Time (days)')

%-- Figures for GFAP growth ----
figure
plot(Vt)
hold on
plot(Vn)
xlabel('Time (days)')
ylabel ('Volume (ng/ml)')
legend('Tumour volume','Necrotic volume')


disp(['current detection time=', num2str(tthresh1)]);
disp(['current detection volume=', num2str(vthresh1)]);
