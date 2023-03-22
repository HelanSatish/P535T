%********************************************************************************************
%               2D TT model propagation using Method of discrete connection of cells
%                             Model Equations of TP06 
%               Helan Satish & M. Ramasubba Reddy - BISP Lab, IIT Madras, India.
%********************************************************************************************


function pdot = Model_Equations(dt,p,nor,noc,iex)

svolt   = p.v;
sm  = p.m;
sh  = p.h;
sj  = p.j;
sd  = p.d;
sf  = p.f;
sf2 = p.f2;
sfcass = p.fcass;
sr = p.r;
ss = p.s;
sxs    = p.xs;
sxr1   = p.xr1;
sxr2  = p.xr2;
sRR  = p.Rprime;
Cai    = p.Cai;
CaSR    = p.Casr;
CaSS    = p.Cass;
Nai    = p.Nai;
Ki    = p.Ki;

%Model Parameters
data=Constants_TP06;
HT=dt;

%Reversal potentials
Ek=data.RTONF*(log((data.Ko./Ki)));
Ena=data.RTONF*(log((data.Nao./Nai)));
Eks=data.RTONF*(log((data.Ko+data.pKNa*data.Nao)./(Ki+data.pKNa.*Nai)));
Eca=0.5*data.RTONF*(log((data.Cao./Cai)));

%Compute currents
%Sodium current
INa=data.GNa*(sm.^3).*sh.*sj.*(svolt-Ena);
AM=1./(1.+exp((-60.-svolt)/5.));
BM=0.1./(1+exp((svolt+35)./5))+0.1./(1+exp((svolt-50)./200));
TAU_M=AM.*BM;
M_INF=1./((1.+exp((-56.86-svolt)./9.03)).*(1.+exp((-56.86-svolt)/9.03)));

sv1=(svolt>-40);
sv2=(svolt<-40);
AH1=zeros(nor,noc);
BH1=(0.77./(0.13.*(1.+exp(-(svolt.*sv1+10.66)./11.1))));BH1=sv1.*BH1;
AJ1=zeros(nor,noc);
BJ1=(0.6.*exp((0.057).*svolt.*sv1)./(1.+exp(-0.1.*(svolt.*sv1+32.))));BJ1=sv1.*BJ1;

AH2=(0.057.*exp(-(svolt.*sv2+80.)/6.8));AH2=sv2.*AH2;
BH2=(2.7.*exp(0.079.*svolt.*sv2)+(3.1e5).*exp(0.3485.*svolt.*sv2));BH2=sv2.*BH2;
AJ2=(((-2.5428e4).*exp(0.2444.*svolt.*sv2)-(6.948e-6).*exp(-0.04391.*svolt.*sv2)).*(svolt.*sv2+37.78)./(1.+exp(0.311.*(svolt.*sv2+79.23))));AJ2=sv2.*AJ2;
BJ2=(0.02424.*exp(-0.01052.*(svolt.*sv2))./(1.+exp(-0.1378.*(svolt.*sv2+40.14))));BJ2=sv2.*BJ2;
AH=AH1+AH2;
BH=BH1+BH2;
AJ=AJ1+AJ2;
BJ=BJ1+BJ2;
TAU_H=1./(AH+BH);
TAU_J= 1./(AJ+BJ);
H_INF=1./((1.+exp((svolt+71.55)./7.43)).*(1.+exp((svolt+71.55)./7.43)));
J_INF=H_INF;

%L-type calcium current
ICaL=data.GCaL.*sd.*(sf.^2).*sfcass.*4.*(svolt-15).*((data.F^2)./(data.R*data.T)).*(0.25.*exp(2.*(svolt-15)*data.F./(data.R*data.T)).*CaSS-data.Cao)./(exp(2.*(svolt-15).*data.F./(data.R*data.T))-1);
D_INF=1./(1.+exp((-8-svolt)./7.5));
Ad=1.4./(1.+exp((-35-svolt)./13))+0.25;
Bd=1.4./(1.+exp((svolt+5)./5));
Cd=1./(1.+exp((50-svolt)./20));
TAU_D=Ad.*Bd+Cd;
F_INF=1./(1.+exp((svolt+20)./7));
Af=1102.5*exp(-(svolt+27).*(svolt+27)./225);
Bf=200./(1+exp((13-svolt)./10.));
Cf=(180./(1+exp((svolt+30)./10)))+20;
TAU_F=Af+Bf+Cf;
F2_INF=0.67./(1.+exp((svolt+35)./7))+0.33;
Af2=600.*exp(-(svolt+25).*(svolt+25)./170);
Bf2=31./(1.+exp((25-svolt)./10));
Cf2=16./(1.+exp((svolt+30)./10));
TAU_F2=Af2+Bf2+Cf2;
FCaSS_INF=0.6./(1+(CaSS./0.05).*(CaSS./0.05))+0.4;
TAU_FCaSS=80./(1+(CaSS./0.05).*(CaSS./0.05))+2.;

%transient outward current
nor=250;
noc=100;
Ito(:,1:25)=data.Gto1.*sr(:,1:25).*ss(:,1:25).*(svolt(:,1:25)-Ek(:,1:25));
Ito(:,26:60)=data.Gto2.*sr(:,26:60).*ss(:,26:60).*(svolt(:,26:60)-Ek(:,26:60));
Ito(:,61:100)=data.Gto3.*sr(:,61:100).*ss(:,61:100).*(svolt(:,61:100)-Ek(:,61:100));
R_INF=1./(1.+exp((20-svolt)/6.));
TAU_R=9.5.*exp(-(svolt+40).*(svolt+40)./1800)+0.8;
S_INF=1./(1.+exp((svolt+20)/5.));
TAU_S=85.*exp(-(svolt+45).*(svolt+45)./320.)+5./(1+exp((svolt-20)./5))+3.;
%%
%slow delayed rectifier current
nor=250;
noc=100;
%Gkss = zeros(nor,noc);
IKs(:,1:25)=data.Gks1.*(sxs(:,1:25).^2).*(svolt(:,1:25)-Eks(:,1:25));
IKs(:,26:60)=data.Gks2.*(sxs(:,26:60).^2).*(svolt(:,26:60)-Eks(:,26:60));
IKs(:,61:100)=data.Gks3.*(sxs(:,61:100).^2).*(svolt(:,61:100)-Eks(:,61:100));
%%
Xs_INF=1./(1.+exp((-5.-svolt)./14));
Axs=(1400./(sqrt(1.+exp((5.-svolt)/6))));
Bxs=(1./(1.+exp((svolt-35.)/15)));
TAU_Xs=Axs.*Bxs+80;

%rapid delayed rectifier current
IKr=data.Gkr*sqrt(data.Ko/5.4).*sxr1.*sxr2.*(svolt-Ek);
Xr1_INF=1./(1.+exp((-26.-svolt)/7.));
axr1=450./(1.+exp((-45.-svolt)/10.));
bxr1=6./(1.+exp((svolt-(-30.))/11.5));
TAU_Xr1=axr1.*bxr1;
Xr2_INF=1./(1.+exp((svolt-(-88.))/24.));
axr2=3./(1.+exp((-60.-svolt)/20.));
bxr2=1.12./(1.+exp((svolt-60.)/20.));
TAU_Xr2=axr2.*bxr2;

%inward rectifier current
Ak1=0.1./(1.+exp(0.06.*(svolt-Ek-200)));
Bk1=(3.*exp(0.0002.*(svolt-Ek+100))+ exp(0.1.*(svolt-Ek-10)))./(1.+exp(-0.5.*(svolt-Ek)));
rec_iK1=Ak1./(Ak1+Bk1);
IK1=data.GK1*sqrt(data.Ko/5.4).*rec_iK1.*(svolt-Ek);

%sodium/calcium exchanger current
INaCa=data.knaca*(1./(data.KmNai^3+data.Nao^3))*(1./(data.KmCa+data.Cao)).* (1./(1+data.ksat.*exp((data.n-1).*svolt.*data.F/(data.R*data.T)))).* (exp(data.n.*svolt.*data.F/(data.R*data.T)).*Nai.^3*data.Cao-exp((data.n-1).*svolt.*data.F/(data.R*data.T))*data.Nao.^3.*Cai.*data.alpha);

%Sodium/potassium pump current
rec_iNaK=(1./(1.+0.1245*exp(-0.1.*svolt*data.F/(data.R*data.T))+0.0353*exp(-svolt.*data.F/(data.R*data.T))));
INaK=data.knak*(data.Ko./(data.Ko+data.KmK)).*(Nai./(Nai+data.KmNa)).*rec_iNaK;

%Plateau potassium current
rec_ipK=1./(1.+exp((25-svolt)./5.98));
IpK=data.GpK*rec_ipK.*(svolt-Ek);

%Plateau Calcium current
IpCa=data.GpCa.*Cai./(data.KpCa+Cai);

%Background currents
IbNa=data.GbNa.*(svolt-Ena);
IbCa=data.GbCa.*(svolt-Eca);

%Calcium dynamics
%update concentrations
kCaSR=data.maxsr-((data.maxsr-data.minsr)./(1+(data.EC./CaSR).*(data.EC./CaSR)));
k1=data.k1prime./kCaSR;
k2=data.k2prime.*kCaSR;
dRR=data.k4.*(1-sRR)-k2.*CaSS.*sRR;
newsRR=sRR+HT.*dRR;
sOO=k1.*CaSS.^2.*sRR./(data.k3+k1.*CaSS.^2);

Irel=data.Vrel*sOO.*(CaSR-CaSS);
Ileak=data.Vleak.*(CaSR-Cai);
Iup=data.Vmaxup./(1.+((data.Kup^2)./(Cai.^2)));
Ixfer=data.Vxfer.*(CaSS-Cai);

%CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
CaCSQN=1./(1+((data.Bufsr*data.Kbufsr)./((CaSR+data.Kbufsr).^2)));
dCaSR=CaCSQN.*(Iup-Irel-Ileak);
newCaSR=CaSR+HT.*dCaSR;

%CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
CaSSBuf=1./(1+((data.Bufss*data.Kbufss)./((CaSS+data.Kbufss).^2)));
dCaSS=CaSSBuf.*(-(Ixfer.*(data.Vc/data.Vss))+(Irel.*(data.Vsr/data.Vss))+(-ICaL.*data.inversevssF2*data.CAPACITANCE));
newCaSS=CaSS+HT.*dCaSS;

%CaBuf=Bufc*Cai/(Cai+Kbufc);
CaBuf=1./(1+((data.Bufc*data.Kbufc)./((Cai+data.Kbufc).^2)));
dCai=CaBuf.*((-(IbCa+IpCa-2*INaCa).*data.inverseVcF2*data.CAPACITANCE)+((Ileak-Iup).*(data.Vsr/data.Vc))+Ixfer);
newCai=Cai+HT.*dCai;

%Sodium/Potassium dynamics
dNai=-(INa+IbNa+3*INaK+3*INaCa)*data.inverseVcF*data.CAPACITANCE;
newNai=Nai+HT.*dNai;

dKi=-(iex+IK1+Ito+IKr+IKs-2*INaK+IpK)*data.inverseVcF*data.CAPACITANCE;
newKi=Ki+HT.*dKi;

%Determine total current
sItot = iex+(IKr+IKs+IK1+Ito+INa+IbNa+ICaL+IbCa+INaK+INaCa+IpCa+IpK);

%Update gates
newsm = M_INF-(M_INF-sm).*exp(-HT./TAU_M);
newsh = H_INF-(H_INF-sh).*exp(-HT./TAU_H);
newsj = J_INF-(J_INF-sj).*exp(-HT./TAU_J);
newsxr1 = Xr1_INF-(Xr1_INF-sxr1).*exp(-HT./TAU_Xr1);
newsxr2 = Xr2_INF-(Xr2_INF-sxr2).*exp(-HT./TAU_Xr2);
newsxs = Xs_INF-(Xs_INF-sxs).*exp(-HT./TAU_Xs);
newss= S_INF-(S_INF-ss).*exp(-HT./TAU_S);
newsr= R_INF-(R_INF-sr).*exp(-HT./TAU_R);
newsd = D_INF-(D_INF-sd).*exp(-HT./TAU_D);
newsf =F_INF-(F_INF-sf).*exp(-HT./TAU_F);
newsf2 =F2_INF-(F2_INF-sf2).*exp(-HT./TAU_F2);
newsfcass =FCaSS_INF-(FCaSS_INF-sfcass).*exp(-HT./TAU_FCaSS);

%update voltage
newsvolt= svolt + dt.*(-sItot);

%pdot=[newsvolt newsm newsh newsj newsd newsf newsf2 newsfcass newsr newss newsxs newsxr1 newsxr2 newsRR newCai newCaSR newCaSS newNai newKi];
pdot=struct('v',newsvolt,'m',newsm,'h',newsh,'j', newsj,'d', newsd,'f',newsf,'f2', newsf2,'fcass', newsfcass,'r', newsr,'s',  newss ,'xs',newsxs,'xr1', newsxr1,'xr2', newsxr2,'Rprime', newsRR,'Cai', newCai,'Casr', newCaSR,'Cass',newCaSS ,'Nai',newNai,'Ki',newKi);

end
