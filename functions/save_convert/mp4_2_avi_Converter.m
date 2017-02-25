%% mp4 2 avi for imageJ
% output path: C:\Users\noambox\Documents\MATLAB

[filename,path] = uigetfile({'*.mp4','MP4 Videos';...
          '*.*','All Files' },'Choose MP4 Video to convert',...
          'D:\# Projects (Noam)\# SLITE\# DATA');
files = dir(fullfile([path,filename]));
for i = 1:length(files)
filename = files(i).name;

    %create objects to read and write the video
    readerObj = VideoReader(files(i).name);
    str = ['AVIconverted_',filename]; % convertion name
    [~,fileBase,fileExt] = fileparts(str);
    writerObj = VideoWriter([fileBase '.avi'],'Uncompressed AVI');
    %open AVI file for writing
    
    open(writerObj);
    %read and write each frame
    for k = 1:readerObj.NumberOfFrames
       img = read(readerObj,k);
       writeVideo(writerObj,img);
    end
    close(writerObj);
end

