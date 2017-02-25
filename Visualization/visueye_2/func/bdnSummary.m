function [rf_3d , rf_3d_blur] = bdnSummary(handles,cell_data)
%bdnSummary(handles,cell_data); returns smoothed spatio-temporal kernel 
% spike 

if(cell_data.props{:}.isnew)
    sigma = 4; % sigma = 4 yields a 17 (=5 pix) filtersize | sigma = 5 yields 21 (~6 pixs)
else
    sigma = 4 * 3; % sigma = 12 yields a 49 (=5 pix) filtersize | sigma = 16 yields 61 (~6 pixs)
end
stim_fs = 5; % hz
rf_bounds = 400; % msec
% if(get(handles.param_SLDR,'Value'))
%    s_thres = 0.01+0.2*get(handles.param_SLDR,'Value');
% else
%     s_thres = 0.01;
% end
if(isfield(handles,'figure1'))
    if(str2double(get(handles.param_ETXT,'String')))
        peri_stim_time = str2double(get(handles.param_ETXT,'String'))+rf_bounds; % in msec
    end
elseif(handles.param_ETXT ~= 0 )
    peri_stim_time = handles.param_ETXT+rf_bounds; % in msec
else
    peri_stim_time = 1300; % in msec
end

pstl = floor(peri_stim_time / (1e3/stim_fs)); % peri-stimulus lagstimulus is in 5 Hz s.t. 1 sample = 200 msec 
l = cell_data.props{1}.signalLenght;
dt = 1/cell_data.props{1}.imaging_rate;
time = linspace(0,(l-1) * dt,l);
if(cell_data.props{:}.isnew)
    stimulus = handles.bdn_src.new;
else
    stimulus = handles.bdn_src.old;
end

% RF mapping
w = cell_data.data.S(cell_data.data.S ~= 0);
if(~any(w))
    return;
end
w = w/max(w); % norma weights
tidx = find(cell_data.data.S ~= 0);
tidx_2 = repmat(tidx,pstl,1);
if(isfield(handles,'offset'))
    offset_samp = handles.offset/ (1e3/stim_fs);
    tidx_3 = [fliplr((1:offset_samp))';-(0:pstl)'];
    tidx_2 = repmat(tidx,pstl+length(1:offset_samp)+1,1);
else
    tidx_3 = -(1:pstl)';
end
tidx = bsxfun(@plus, tidx_2, tidx_3);
tidx_remove = find(tidx(1,:) > l | tidx(end,:) < 1);
tidx(:, tidx_remove) = [];
w(:, tidx_remove) = [];

% map tidx to stimulus samples
resampling_ratio = cell_data.props{:}.imaging_rate/stim_fs;
rs_vec = repmat((1:ceil(l/resampling_ratio )),2,1);
rs_vec = rs_vec(:);
rs_vec = rs_vec(1:l);

% build the matrix container
% disp('DEBUG MODE');
% tidx = randi(length(cell_data.data.S), pstl, 300); % FIXME test for random picks of spikes

mat_container = zeros(size(stimulus.M,1),size(stimulus.M,2),size(tidx,1),length(tidx));
for k = 1:size(tidx,2)
    try
        mat_container(:,:,:,k) =  2 .* w(k) .* stimulus.M(:,:,rs_vec(tidx(:,k)));
    catch exception
        disp((tidx(:,k)));
    end
end
rf_3d = mean(mat_container,4);
d_samp = (pstl-rf_bounds/ (1e3/stim_fs));
rf_3d_blur = zeros(size(mat_container,1),size(mat_container,2),size(tidx,1));
for k = 1:(size(tidx,1))
    rf_3d_blur(:,:,k) = imgaussfilt(rf_3d(:,:,k),sigma);
end

if(~(isfield(handles,'figure1')))
    return;
end

% Plots
% TODO: take this out the function

map = [[linspace(0,0,24)';linspace(0,1,26)'] [linspace(0,0,40)';linspace(0,1,10)'] linspace(0,0.6,50)'];
axes(handles.axes2);
colormap(map);
mat_now = std(rf_3d_blur(:,:,1:d_samp),0,3);
him = imagesc(mat_now);
ax_position = fliplr(size(him.CData));
[x_cell, y_cell] = findCircleCenter(cell_data.data.location(1), cell_data.data.location(2), round(ax_position),cell_data.props{:}.isnew);
if(cell_data.props{:}.isnew)
    circles = int32([x_cell y_cell 1;x_cell y_cell 2]);
else
    circles = int32([x_cell y_cell 1;x_cell y_cell 4]);
end

% [X,Y] = meshgrid((0:250),(0:150));
% T2 = [[0.18003081420734915, -0.0037810874601845275, 0]' [-8.197621367610952e-4, -0.18047912349874015, 0]' [54.859700082945054, 100.85998914002913,1]'];
% [s,k] = deal(ax_position(1)/200, ax_position(2)/150);
% coord_t = (T2 *[X(:), Y(:), ones(numel(X),1)]'); 
% X1 = round( coord_t(1,:) * s ); 
% Y1 = round( coord_t(2,:) * k );
% circles = int32([X1' Y1' ones(numel(X1),1)]);
% figure();
% colormap(map);
% J = step(si , mat_now, circles);
% imagesc(J)

si = updateShapeInserter(max(mat_now(:)));
J = step(si, std(rf_3d_blur(:,:,1:d_samp),0,3), circles);
set(gca,'NextPlot','replacechildren');
set(gca, 'Ydir', 'reverse');
set(gca, 'YAxisLocation', 'Right');
if(cell_data.props{:}.isnew);J = imcrop(J,[144/4 216/4-4 348/4 (264)/4+8]);end;
imagesc(addScaleBar(J,cell_data.props{:}.isnew, max(mat_now(:))));
title(['STA latency (msec) - ',num2str(peri_stim_time - rf_bounds)]);
axis tight;

figure(1); clf;
colormap(map);
% colormap('cool');colormap(flipud(colormap));
set(gcf,'position',[2000 100 1000 800]);
subplot(2,2,1)
mat_now = std(rf_3d_blur(:,:,(1:d_samp-2)),0,3);
si = updateShapeInserter(max(mat_now(:)));
J = step(si , mat_now, circles);
if(cell_data.props{:}.isnew);J = imcrop(J,[144/4 216/4-4 348/4 (264)/4+8]);end;
imagesc(addScaleBar(J,cell_data.props{:}.isnew,max(mat_now(:))));
title(['STA latency (msec) - ',num2str(peri_stim_time - rf_bounds - 400)]);

subplot(2,2,2)
mat_now = std(rf_3d_blur(:,:,(1:d_samp-1)),0,3);
si = updateShapeInserter(max(mat_now(:)));
J = step(si , mat_now, circles);
if(cell_data.props{:}.isnew);J = imcrop(J,[144/4 216/4-4 348/4 (264)/4+8]);end;
imagesc(addScaleBar(J,cell_data.props{:}.isnew,max(mat_now(:))));
title(['STA latency (msec) - ',num2str(peri_stim_time - rf_bounds - 200)]);

subplot(2,2,3)
mat_now = std(rf_3d_blur(:,:,(1:d_samp+1)),0,3);
si = updateShapeInserter(max(mat_now(:)));
J = step(si , mat_now, circles);
if(cell_data.props{:}.isnew);J = imcrop(J,[144/4 216/4-4 348/4 (264)/4+8]);end;
imagesc(addScaleBar(J,cell_data.props{:}.isnew,max(mat_now(:))))
title(['STA latency (msec) - ',num2str(peri_stim_time - rf_bounds + 200)]);

subplot(2,2,4)
mat_now = std(rf_3d_blur(:,:,(1:d_samp+2)),0,3);
si = updateShapeInserter(max(mat_now(:)));
J = step(si , mat_now, circles);
if(cell_data.props{:}.isnew);J = imcrop(J,[144/4 216/4-4 348/4 (264)/4+8]);end;
imagesc(addScaleBar(J,cell_data.props{:}.isnew,max(mat_now(:))))
title(['STA latency (msec) - ',num2str(peri_stim_time - rf_bounds + 400)]);

end

function [x,y] = findCircleCenter(x0,y0, h_position, isnew)
% translate cell location according to estimate form OneNote

T1 = [[0.8956347255271192, -0.006083132214985391, 0]' [0.005243919072709428, -0.662675926149547, 0]' [160.1354410707739, 421.5876157241306, 1]']; % custom made
if(isnew)
%     T2 = [[0.5426369331949202, 6.417199507468235e-5, 0]' [5.326217994102055e-6, 0.5710082471585562, 0 ]' [72.71713363542452, 30.962022148317914, 1]'];
    T2 = [[0.18003081420734915, -0.0037810874601845275, 0]' [-8.197621367610952e-4, -0.18047912349874015, 0]' [54.859700082945054, 100.85998914002913,1]'];
    [w,h] = deal(h_position(1)/200, h_position(2)/150);
else
%     T2 = [[1.247079475972392, 0.0029743930210039847,0]' [-0.0029743930210039847, 1.247079475972392,0]' [-101.82311986332549, -227.32686356053898,1]'];
    T2 = [[0.4052384483847369, 0.004615168124209535,0]' [-0.012652302229976614, -0.39706379697304445,0]' [47.304568570397166, 105.67567123777455,1 ]'];
    [w,h] = deal(h_position(1)/200, h_position(2)/150);
end

% coord_t = (T2 * (T1 * [x0, y0, 1]')); 
coord_t = (T2 *[x0, y0, 1]'); 
x = round( coord_t(1) * w ); 
y = round( coord_t(2) * h );
end

function [si] = updateShapeInserter(val)
    si = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',val);
end

function [im] = addScaleBar(im, isnew, val)
if(isnew)
    im((70:73),(3:23)) = val; % 45 pix are 490 um -> 2*9.18 pix for 200 um
else
    im((140:145),(5:45)) = val; % 100 pix are 490 um -> 2*20.408 pix for 200 um
end
end
