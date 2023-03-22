%********************************************************************************************
%               2D TT model propagation using Method of discrete connection of cells
%               Extract frames from the video
%               Helan Satish & M. Ramasubba Reddy - BISP Lab, IIT Madras, India.
%********************************************************************************************
%% frames
ReadObj = VideoReader('P535T_HM.mp4'); 
CurFrame = 0;
GetFrame = [6401 6862 7397 7468 7504 7642 7858 7985 8163 8347 8379];
while hasFrame(ReadObj)
    CurImage = readFrame(ReadObj);
    CurFrame = CurFrame+1;
    if ismember(CurFrame, GetFrame)
        imwrite(CurImage, sprintf('frame%d.tif', CurFrame));
    end
end