%********************************************************************************************
%               2D TT model propagation using Method of discrete connection of cells
%               Saving the frames as video
%               Helan Satish & M. Ramasubba Reddy - BISP Lab, IIT Madras, India.
%********************************************************************************************
function [] = video(vfin)
dt=0.05;
[a1 b1 c1]= size(vfin);
vidObj = VideoWriter('tmvideo1.avi');
     open(vidObj);
        for k= 1:1:c1
            k
            xx=vfinn(:,:,k);
            ih=imagesc(xx);caxis([-85 50]);colormap(hot);colorbar;
            hold on
            line([25 25 ], [0 250 ],'color','w','LineWidth',2)
            line([60 60], [0 250 ],'color','w','LineWidth',2)
            hold on
            xlabel('Number of Rows');
            ylabel('Number of Columns')
            title(sprintf('Time : %i milliseconds ',round(k*dt)))
            mov=getframe;
            k=k+1;
            writeVideo(vidObj, mov);
            hold on
 end       
close(vidObj);
end