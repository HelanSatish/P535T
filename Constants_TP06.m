%********************************************************************************************
%               2D TT model propagation using Method of discrete connection of cells
%               Constants for TP06 2d
%               Helan Satish & M. Ramasubba Reddy - BISP Lab, IIT Madras, India.
%********************************************************************************************
function data = Constants_TP06()

% General parameters
data.F=96485;   % Faraday constant
data.R=8314;    % Gas constant
data.T=310;     % Temperature
data.z=1;data.zca=2;
data.CAPACITANCE=0.185;
data.Vc=0.016404;
data.Vsr=0.001094;
data.Vss=0.00005468;
data.RTONF=(data.R*data.T)/data.F;
data.inverseVcF2=1/(2*data.Vc*data.F);
data.inverseVcF=1./(data.Vc*data.F);
data.inversevssF2=1/(2*data.Vss*data.F);

% Sodium
data.Nao=140;
data.GNa=14.838;

% Calcium
data.Cao=2;

% Potassium
data.Ko=5.4;
data.GK1=5.405;
%% Gto and Gks
data.Gto1=0.073;
data.Gto2=0.294;
data.Gto3=0.294;
data.Gkr=0.153;
% P535T: 
%WT: 0.392 0.098 0.392
%HT: 0.242 0.060 0.242
%HM: 0.159 0.039 0.159
%% WT
data.Gks1=0.392;%epi/endo
data.Gks2=0.098;
data.Gks3=0.392;
%% HT
% data.Gks1=0.242;
% data.Gks2=0.060;
% data.Gks3=0.242;
%% HM
% data.Gks1=0.159;
% data.Gks2=0.039;
% data.Gks3=0.159;
%%
% Other parameters
data.pKNa=0.03;
data.GCaL=3.98e-5;
data.knaca=1000;
data.n=0.35;
data.KmCa=1.38;
data.KmNai=87.5;
data.ksat=0.1;
data.alpha=2.5;
data.knak=2.724;
data.KmK=1;data.KmNa=40;
data.GpCa=0.1238;data.GpK=0.0146;data.KpCa=0.0005;
data.GbNa=0.00029;data.GbCa=0.000592;
data.Vmaxup=0.006375;
data.Kup=0.00025;
data.Vrel=0.102;
data.k1prime=0.15;
data.k2prime=0.045;
data.k3=0.06;
data.k4=0.005;
data.EC=1.5;
data.minsr=1;data.maxsr=2.5;
data.Vleak=0.00036;
data.Vxfer=0.0038;
data.Bufc=0.2;data.Kbufc=0.001;
data.Bufsr=10;data.Kbufsr=0.3;
data.Bufss=0.4;data.Kbufss=0.00025;
end

