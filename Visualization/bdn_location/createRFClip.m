function [] = createRFClip(fname, start_frame, nFrames, rfm_blurr, time, map)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

fig = figure(6);
h1 = axes;
vw = VideoWriter(fname,'Motion JPEG AVI');
vw.Quality = 60;
vw.FrameRate = 1;
vw.open();
for k = fliplr(start_frame:nFrames)
    
    % plot stimulus frame
    axes(h1);colormap(map);imagesc(rfm_blurr(:,:,k));axis(h1,'tight');
    title(h1,['Time before stimulus - ',num2str(-time(k)),' [msec]'])
    set(gca,'FontSize',16);

    frame = getframe(gcf);
    writeVideo(vw,frame);
    drawnow;
end
vw.close();
close(fig);
end

