%********************************************************************************************
%               2D TT model propagation using Method of discrete connection of cells
%                  Helan Satish & M. Ramasubba Reddy - BISP Lab, IIT Madras, India.
%********************************************************************************************

function [iinew] = Potential(v)
%% ######################################
    %Update v using 8 neighbours
    %Horizontal & Vertical components
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
    
    %Total current
    iinew=((vxx+vd1+vd2).*Gxx)+(vyy.*Gyy);
end