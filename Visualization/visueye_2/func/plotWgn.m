function [ ] = plotWgn( handles,time,kernel ,all_signals, stimulus, cell_data )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
cla(handles.axes2);
axes(handles.axes2);
set(handles.axes2,'Visible','On');
set(handles.axes2,'NextPlot','replacechildren');
for k = 1:size(all_signals,2)
    all_signals(:,k) = mapminmax_nc(all_signals(:,k),max(kernel),min(kernel));
end
% plot(time,all_signals,'Color',[0.85    0.85    0.85]);
hold on;
plot(time, kernel,'.--','LineWidth',1);
xlabel('time [sec]');ylabel('STA intensity');
smooth_kernel = smooth(kernel, 0.08);
if(any(kernel))
    plot(time, smooth_kernel,'LineWidth',1.5);
    axis([time(end) time(1) (min(kernel) - std(kernel)) (max(kernel) + std(kernel))])
    line([0 0], [(min(kernel) - 2*std(kernel)) (max(kernel) + 2*std(kernel))],'LineStyle','--','Color',[0.4941    0.4941    0.4941]);
end
hold off;
set(gca,'FontSize',11);

fig1 = figure(1);clf;
set(fig1,'position',[ 2000 400 600 600]);
ax1 = subplot(10,6,[(1:11),(21:30)]);
copyobj(allchild(handles.axes2),ax1);
xlabel('time [sec]');ylabel('STA intensity');
axis([time(end) time(1) (min(kernel) - std(kernel)) (max(kernel) + std(kernel))])

ax2 = subplot(10,6,[(41:10:50),(51:56)]);
cla(ax2);
copyobj(allchild(handles.axes1),ax2);
xlabel('Time [sec]');ylabel('norma Ca^{2+} est.');
axis([get(handles.axes1,'XLim') get(handles.axes1,'YLim')]);

%% plot simulated linear rate process comparison
% convolve estimated kernel with the outcome of previus step
kernel_trunc = smooth_kernel(time'<=0);
dt_stim = 1/30; % 30 Hz
dt_r = 1/cell_data.props{:}.imaging_rate;
time_stimulus = (0:dt_stim:(length(kernel_trunc)-1)*dt_r);
time_r = (0:dt_r:(length(kernel_trunc)-1)*dt_r);

%resampling kernel
kernel_intrp = interp1(time_r, kernel_trunc, time_stimulus,'pchip');
L_step = conv(stimulus, (kernel_intrp),'same'); % already fliiped kernel
lambda_rate = L_step;

% non_linear and truncating
NL_step = L_step;
% NL_step = NL_step.^2;
lambda_rate = NL_step;
lambda_rate = lambda_rate - mean(lambda_rate);
lambda_rate(lambda_rate<=0) = eps;
nl = 0;

% compare to the rate process as an Linear model
figure(6);
time_s = (0:dt_stim:(length(stimulus)-1)*dt_stim); % in sec
l_r = cell_data.props{:}.signalLenght;
time_r = (0:dt_r:(l_r-1)*dt_r);
% L_intrp_blurr = interp1(time_s, L_step_blurr, time_r);
L_intrp = interp1(time_s, lambda_rate, time_r);
response = cell_data.data.S;
% plot(time_r,L_intrp-mean(L_intrp));hold on;plot(time_r,mapminmax_nc(response,0,max(L_intrp(100:end)-mean(L_intrp(100:end)))));hold off;
plot(time_r,L_intrp);hold on;plot(time_r,mapminmax_nc(response,eps,max(L_intrp(100:end)-mean(L_intrp(100:end)))));hold off;
xlabel('Times [sec]');
legend('estimated rate with model','spikes');
if(~nl)
    title('Linear model rate estimation vs. estimated spike process')
else
    title('LN model rate estimation vs. estimated spike process')
end
axis([time_r(1) time_r(end) -inf inf])
set(gca,'FontSize',14);
set(gcf,'Position',[1 100 1200 250]);

% folding squaring and rectifying
cell_data.data.S = cell_data.data.S
figure(7);
par = cell_data.stim{:}.partition;
c_fold = cell_data.data.C;

if length(cell_data.data.C) == length(cell_data.data.S)
    s_fold = cell_data.data.S;
else
    s_fold = cell_data.data.full_S;
end
s_fold = arrayfun(@(x) s_fold(par(x,1):par(x,2)),(1:size(par,2)),'UniformOutput',false);
c_fold = arrayfun(@(x) c_fold(par(x,1):par(x,2)),(1:size(par,2)),'UniformOutput',false);
subplot(211);plot(c_fold{1}); hold on; plot(c_fold{2});hold off; legend('1','2');
subplot(212);plot(s_fold{1});hold on; plot(s_fold{2});hold off;legend('1','2');
end

