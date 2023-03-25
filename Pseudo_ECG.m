%********************************************************************************************
%                                 Generation of Pseudo ECG
%               Helan Satish & M. Ramasubba Reddy - BISP Lab, IIT Madras, India.
%********************************************************************************************

function [L2_norm,times] = Pseudo_ECG(vfin)
nor=250;noc=100;
A=size(vfin);
leng=A(3);%length(vfin);%160;%
vm1=zeros(nor*noc,leng);
r3=1;
for r1=1:nor
    for r2=1:noc
        vm1(r3,:)=reshape(vfin(r1,r2,:),1,leng);        
        r3=r3+1;
    end
end
%%
[row,col,t]=size(vfin);
magR1(1:row,1:col)=zeros;
magR2(1:row,1:col)=zeros;
magR3(1:row,1:col)=zeros;
for i=1:row
    for j=1:col
        D=[i, j, 0];
        
        %Contribution of each dipole to V1
        V1R=[50, 75, 50]; % 25, 95, 50
        R1=V1R-D;
        magR1(i,j)=norm(R1);
                
        %Contribution of each dipole to V2
        V2R=[75,10,50]; % 75, 10, 50
        R2=V2R-D;
        magR2(i,j)=norm(R2);
                
        %Contribution of each dipole to V3
        V3R=[125,80,50]; % 225,100 50
        R3=V3R-D;
        magR3(i,j)=norm(R3);
       
    end
end

for tt=1:t;
    v=vfin(:,:,tt);
    vv=[[v(1,1) v(1,:) v(1,end)];[v(:,1) v v(:,end)];[v(end,1) v(end,:) v(end,end)]];
    vxx=(vv(2:end-1,1:end-2) + vv(2:end-1,3:end) -2*v);
    vyy=(vv(1:end-2,2:end-1) + vv(3:end,2:end-1) -2*v);
    %Diagonal1
    c1=[v(2:end,1);0]; c2=[0;v(1:end-1,end)];
    c=[c1 v c2];
    vv1=[[v(1,:) 0 0];c;[0 0 v(end,:)]];
    vd1=(vv1(1:end-2,1:end-2) + vv1(3:end,3:end) -2*v);
    %Diagonal2
    c3=[0; v(1:end-1,1)]; c4=[v(2:end,end);0];
    c5=[c3 v c4];
    vv2=[[0 0 v(1,:)];c5;[v(end,:) 0 0]];
    vd2=(vv2(3:end,1:end-2) + vv2(1:end-2,3:end) -2*v);
    
    V11mat=vxx./(magR1.^2);%horizontal contribution
    V12mat=vyy./(magR1.^2);%vertial
    V13mat=vd1./(magR1.^2);%downward diagonal
    V14mat=vd2./(magR1.^2);%upward diagonal
    V1(tt)=sum(sum(V11mat+V12mat+V13mat+V14mat));
    
    
    V21mat=vxx./(magR2.^2);%horizontal contribution
    V22mat=vyy./((magR2.^2));%vertial
    V23mat=vd1./((magR2.^2));%downward diagonal
    V24mat=vd2./((magR2.^2));%upward diagonal
    V2(tt)=sum(sum(V21mat+V22mat+V23mat+V24mat));
    
    
    V31mat=vxx./((magR3.^2));%horizontal contribution
    V32mat=vyy./((magR3.^2));%vertial
    V33mat=vd1./((magR3.^2));%downward diagonal
    V34mat=vd2./((magR3.^2));%upward diagonal
    V3(tt)=sum(sum(V31mat+V32mat+V33mat+V34mat));
end
L1=V2-V1;
L2=V3-V1;
L3=V3-V2;
L1_norm=L1./max(abs(L2));
L2_norm=L2./max(abs(L2));
L3_norm=L3./max(abs(L2));

%% Plot
time=1:leng;
timest=time*100*0.05/1000;
figure();
ax = axes;
grid(ax);
grid(ax, 'minor');
ax.GridColor = [ 1 0 0 ];
ax.MinorGridColor = [1 0 0];
hold on;
plot(timest,L2_norm(:,1:leng),'k','linewidth',2.5);
set(gca,'FontName','Times New Roman','FontSize',20);
box on;
xlabel('Time (s)');
ylabel ('V_m');
end
