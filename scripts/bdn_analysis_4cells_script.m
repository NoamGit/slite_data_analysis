% investiagte BDN
dt_stim = (1e3/5);
B = load('C:\Users\noambox\Documents\Sync\Neural data\bdn_focus.mat'); 
bdn_focus = B.bdn_focus;
stim_src_new = load('C:\Users\noambox\Documents\Sync\stimulus files\src\BDN_stimulus_new_5hz.mat'); 
stim_src_old = load('C:\Users\noambox\Documents\Sync\stimulus files\src\BDN_stimulus_5hz.mat');
param.bdn_src.new = stim_src_new;
param.bdn_src.old = stim_src_old;
scale_old = 200/24;
y_offset = 25;
stim = param.bdn_src.old; % all 4 examples are old
    
% cell specifics
    %% C64
    ii = 1; 
    param.param_ETXT = 1200;
    param.offset = 200;
    pstl = floor(param.param_ETXT / dt_stim); % peri-stimulus lagstimulus is in 5 Hz s.t. 1 sample = 200 msec 
    offstl = (param.offset/dt_stim);
    f_map_start = offstl+2;
    d_samp = (pstl - 200/dt_stim)-2;
    fname = 'C64.avi'; start_frame = 2; 
    %% C15
    ii = 2; 
    param.param_ETXT = 1200;
    param.offset = 600;
    pstl = floor(param.param_ETXT / dt_stim); % peri-stimulus lagstimulus is in 5 Hz s.t. 1 sample = 200 msec 
    offstl = (param.offset/dt_stim);
    f_map_start = offstl+2;
    d_samp = (pstl - 200/dt_stim)-2;
    fname = 'C15.avi'; start_frame = 2;  
    %% C47
    ii = 3;
    fname = 'C47.avi'; start_frame = 2; 
%     ii = 5;
%     fname = 'C47_2.avi'; start_frame = 2; 
    param.param_ETXT = 800;
    param.offset = 400;
    pstl = floor(param.param_ETXT / dt_stim); % peri-stimulus lagstimulus is in 5 Hz s.t. 1 sample = 200 msec 
    offstl = (param.offset/dt_stim);
    f_map_start = offstl+2;
    d_samp = (pstl -200/dt_stim)-1;
    %% C54
    ii = 4; 
    param.param_ETXT = 1000;
    param.offset = 400;
    pstl = floor(param.param_ETXT / dt_stim); % peri-stimulus lagstimulus is in 5 Hz s.t. 1 sample = 200 msec 
    offstl = (param.offset/dt_stim);
    f_map_start = offstl+2;
    d_samp = (pstl - 200/dt_stim)+1;
    fname = 'C54.avi'; start_frame = 2;  
%% extract RF data for cell 

[x_cell, y_cell] = findCircleCenter(bdn_focus(ii,:).data.location(1), bdn_focus(ii,:).data.location(2),...
    round([200 150]),bdn_focus(ii,:).props{:}.isnew);
map = [[linspace(0,0,24)';linspace(0,1,26)'] [linspace(0,0,40)';linspace(0,1,10)'] linspace(0,0.6,50)'];
[rfm, rfm_blurr] = bdnSummary(param, bdn_focus(ii,:));
mat_now = std(rfm_blurr(:,:,(f_map_start:f_map_start+(d_samp))),0,3);
figure(2);imagesc(mat_now);
title(['STA latency ',num2str(d_samp.*dt_stim),' [msec]']);
colormap(map);
%% RF - extract kernel

params.showflag = true;
params.gfit = [];
[kernel, gfit] =  findKernelFromIm(mat_now, scale_old, [12 24], rfm_blurr, params);
[xi,yi] = meshgrid(1:200,1:150);
clear opts;
opts.tilted = true;
results = autoGaussianSurf(xi,yi,mat_now,opts);
gfit.angle = results.angle;
subplot(1,2,1);imagesc(mat_now);
subplot(1,2,2);imagesc(xi(:),yi(:),results.G);


% plot kernel
figure(4);clf;
set(gcf,'Position',[2016         650         560         342])
kernel = kernel-mean(kernel);
time = -(-offstl * dt_stim:dt_stim:(size(rfm_blurr,3) - offstl-1)*dt_stim);
hold on;
plot(fliplr(time), fliplr(kernel'),'LineWidth',1.5);
xlabel('Time [msec]');ylabel('STA intensity');
axis([time(end) time(1) (min(kernel) - 0.5*std(kernel)) (max(kernel) + 0.5*std(kernel))])
line([0 0], [(min(kernel) - 2*std(kernel)) (max(kernel) + 2*std(kernel))],'LineStyle','--','Color',[0.4941    0.4941    0.4941]);
hold off;
title('RF centered temporal kernel');
set(gca,'FontSize',10);
%% Simulate rate process

% create blurred stimuli
stim_blurr = zeros(size(stim.M));
sigma = 4*3;
nFrames = size(stim.M,3);
for k = 1:nFrames
    stim_blurr(:,:,k) = imgaussfilt(stim.M(:,:,k),sigma);
end

% use estimated RF mask and dotproduct it with the stimulus
mask = mapminmax_nc(imresize(gfit.G, scale_old), 0, 1);
masked_stim = sum(bsxfun(@times,double(reshape(stim.M(y_offset+1:end-y_offset,:,:),numel(mask),[])),reshape(mask,[],1)),1);
masked_blurr_stim = sum(bsxfun(@times,double(reshape(stim_blurr(y_offset+1:end-y_offset,:,:),numel(mask),[])),reshape(mask,[],1)),1);

% convolve estimated kernel with the outcome of previus step
kernel_trunc = kernel(time'<=0);
L_step = conv(masked_stim, (kernel_trunc),'same'); % already fliiped kernel
L_step_blurr = conv(masked_blurr_stim, mapminmax_nc(kernel_trunc,0,max(masked_blurr_stim)),'same');

% compare to the rate process as an Linear model
figure(6);
time_s = (0:dt_stim/1e3:(nFrames-1)*dt_stim/1e3); % in sec
dt_r = 1/bdn_focus(ii,:).props{:}.imaging_rate;
l_r = bdn_focus(ii,:).props{:}.signalLenght;
time_r = (0:dt_r:(l_r-1)*dt_r);
L_intrp_blurr = interp1(time_s, L_step_blurr, time_r);
L_intrp = interp1(time_s, L_step, time_r);
response = bdn_focus(ii,:).data.S;
response(response<0.01) = nan;
plot(time_r,L_intrp_blurr-mean(L_intrp_blurr));hold on;stem(time_r,mapminmax_nc(response,0,max(L_intrp_blurr(100:end)-mean(L_intrp_blurr(100:end)))),'ShowBaseLine','Off' );hold off;
xlabel('times [sec]');
title('linear model rate estimation vs. estimated spike process')
%% ROI - centered window dynamics

% use moved gaussian profile
xy_scaled = [x_cell, y_cell - y_offset]./scale_old;
gfit_drift = gfit;
gfit_drift.x0 = xy_scaled(1);
gfit_drift.y0 = xy_scaled(2);
params.gfit = gfit_drift;
[kernel_window] = findKernelFromIm(mat_now, scale_old, [12 24], rfm_blurr, params);
params.gfit = [];

% plot kernel
figure(5);clf;
set(gcf,'Position',[2016         200         560         342])
kernel_window = kernel_window-mean(kernel_window);
time = -(-offstl * dt_stim:dt_stim:(size(rfm_blurr,3) - offstl-1)*dt_stim);
hold on;
plot(fliplr(time), fliplr(kernel_window'),'LineWidth',1.5);
xlabel('Time [msec]');ylabel('STA intensity');
axis([time(end) time(1) (min(kernel) - 0.5*std(kernel)) (max(kernel) + 0.5*std(kernel))])
line([0 0], [(min(kernel) - 2*std(kernel)) (max(kernel) + 2*std(kernel))],'LineStyle','--','Color',[0.4941    0.4941    0.4941]);
hold off;
title('Cell centered temporal kernel');
set(gca,'FontSize',10);
%% RF - clips

path = 'C:\Users\noambox\Documents\Sync\code\Data Analysis\Visualization\bdn_location';
cd([path]);
nFrames = size(rfm_blurr,3);
createRFClip(fname, start_frame,nFrames, rfm_blurr, time, map);
%% Draft code

% % find gaussian fit 1
% imsize = [150, 200];
% fun = @(p) sum(norm_nc(gauss2d(p(1), p(2), p(3), imsize),6) - norm_nc(mat_now,6)).^2;
% [~,I] = max(mat_now(:));
% [x0,y0] = ind2sub(size(mat_now),I);
% [d_s,d_mx,d_my] = meshgrid((12:30)',(47:105),(48:147));
% gauss_p = patternsearch(fun,[15, x0, y0],[],[],[],[],[10,47,48],[30,105,147]);
% figure(6);imagesc(gauss2d(gauss_p(1),gauss_p(2),gauss_p(3),imsize));
% figure(7);imagesc(mat_now);

