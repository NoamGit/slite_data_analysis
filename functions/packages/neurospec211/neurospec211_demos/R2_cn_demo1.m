% Script: R2_cn_demo1.m
% Script to demonstrate use of NeuroSpec 2.1
% R2 directionality analysis on spike train data from simulated Cortical Neuron networks.
%
% Copyright 2015, David M. Halliday.
% This file is part of NeuroSpec.
%
%    NeuroSpec is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    NeuroSpec is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with NeuroSpec; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%    NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
%

% This script uses the simulated cortical neuron spike train data reported in JIN article:
%  Halliday DM (2015) Nonparametric directionality measures for time series and point process data,
%  Journal of Integrative Neuroscience, 14(2), In Press. DOI: 10.1142/S0219635215300127.

% Load data
load JIN_dat1
% Contents of MAT file:
%  cn01_3a, cn01_3b, cn01_3c, cn01_3d, cn01_3e, cn01_3f, cn01_3a_d, cn01_3e_d:
%   Data matrices for configuration a), ..., h), respectively. Matrices *_d have additional delays.
%  filstem_list   Text identifier for each configuration
%        n_runs   No of runs for each configuration (10)
%          rate   Sampling rate for data (1000 samples/sec)
%           sec   No of seconds of data, data length is (rate*sec)
%     sp_Pn_all   Matrix of Pn values (mean rate) for each point process spike train
%    sp_tot_all   Matrix of spike totals for each point process spike train
%
%  Data matrices are size [rate*sec,3,n_runs] with spike trains stored in 0/1 format

% Set segment length for analysis, using power of 2
% T=2^10 = 1024. 1024 points @ 1000/sec sampling gives segment length of 1.024 sec
% Frequency resolution is inverse of this: 0.977 Hz.
seg_pwr=10;
seg_size=2^seg_pwr;

%-----------------------------------------------------------------------------
% Storage for results from neuron pairs 1-2 and 1-3
n_set=length(filstem_list);
% R2 frequency domain matrices - 12 columns returned from sp2a2_R2
f12=zeros(seg_size/2,12,n_runs,n_set);
f13=zeros(seg_size/2,12,n_runs,n_set);
% R2 time domain matrices - 3 columns returned from sp2a2_R2
t12=zeros(seg_size,3,n_runs,n_set);
t13=zeros(seg_size,3,n_runs,n_set);

% R2 analysis of pair-wise interactions in each set
for set_no=1:n_set
  fil_stem=filstem_list{set_no};
  % select single set of data to use
  dat=eval(fil_stem);
  disp(['File no: ',num2str(set_no),' - ',fil_stem])

  % R2 analysis for individual runs
  for ind_run=1:n_runs
    % Neurons 1-2
    [f12(:,:,ind_run,set_no),t12(:,:,ind_run,set_no),cl12(ind_run,set_no)]=sp2a2_R2(dat(:,1,ind_run),dat(:,2,ind_run),rate,seg_pwr);
    cl12(ind_run,set_no).what=[fil_stem,': 1-2, run: ',num2str(ind_run)];
    cl12(ind_run,set_no).N1=sp_tot_all(ind_run,1,set_no);  % Append spike counts
    cl12(ind_run,set_no).N2=sp_tot_all(ind_run,2,set_no);
    cl12(ind_run,set_no).P1=sp_Pn_all(ind_run,1,set_no);   % Append expected values for Point Process spectra
    cl12(ind_run,set_no).P2=sp_Pn_all(ind_run,2,set_no);
    
    % Neurons 1-3
    [f13(:,:,ind_run,set_no),t13(:,:,ind_run,set_no),cl13(ind_run,set_no)]=sp2a2_R2(dat(:,1,ind_run),dat(:,3,ind_run),rate,seg_pwr);
    cl13(ind_run,set_no).what=[fil_stem,': 1-3, run: ',num2str(ind_run)];
    cl13(ind_run,set_no).N1=sp_tot_all(ind_run,1,set_no);  % Append spike counts
    cl13(ind_run,set_no).N2=sp_tot_all(ind_run,3,set_no);
    cl13(ind_run,set_no).P1=sp_Pn_all(ind_run,1,set_no);   % Append expected values for Point Process spectra
    cl13(ind_run,set_no).P2=sp_Pn_all(ind_run,3,set_no);
  end 
end

%-----------------------------------------------------------------------------
% Generate table 1 in JIN
max_freq=250;  % Metrics over restricted frequency range, this f_alpha value

% Storage for R2 metrics
table_dat12=zeros(8,4,10);
table_dat13=zeros(8,4,10);
for set_no=1:n_set
  for ind_run=1:n_runs
    % Calculate R2 metrics over restricted frequency range
    [cl12a_f(ind_run,set_no),table_dat12(set_no,:,ind_run)]=R2_f1(f12(:,:,ind_run,set_no),cl12(ind_run,set_no),max_freq);
    [cl13a_f(ind_run,set_no),table_dat13(set_no,:,ind_run)]=R2_f1(f13(:,:,ind_run,set_no),cl13(ind_run,set_no),max_freq);
  end 
end

% Display mean values over 10 runs
disp(' ')
disp('mean R2 for 1->2:')
disp('Overall, Zero Lag, Forward,  Reverse')
disp(num2str(mean(table_dat12,3),3))
disp(' ')
disp('mean R2 for 1->3:')
disp('Overall, Zero Lag, Forward,  Reverse')
disp(num2str(mean(table_dat13,3),3))

% Generate Table 1 for JIN manuscript, uses means over 10 runs and percentage of R2 in each direction
R2_table=zeros(8,10);
R2_table(:,[1 2 4 6 7 9])=[mean(table_dat12(:,[1 3 4 ],:),3) mean(table_dat13(:,[1 3 4 ],:),3)];

% Add percentages for forward and reverse metrics
R2_table(:,[3 5 8 10])=round(100*R2_table(:,[2 4 7 9])./R2_table(:,[1 1 6 6] ));

% Display to match Table 1 in JIN
disp(' ')
disp('R^2 table (Table 1 in JIN): ')
disp(num2str(R2_table,'%8.4f %8.4f %4d %8.4f %4d %8.4f %8.4f %4d %8.4f %4d'))

%-----------------------------------------------------------------------------
% Plotting for neurons 1->2, frequency domain, Figure 2 in JIN
freq_all=[250 250 250 100 100 250 250 100];
ch_max_all=[0.3 0.3;0.3 0.25;0.3 0.3;0.15 0.3; 0.15 0.6;0.3 0.3;0.3 0.3;0.150 0.6];
label_text=['a','b','c','d','e','f','g','h'];

% Landscape layout
freq_label_ind=4;
line_width=0.5;
subplot_1=2;
subplot_2=4;

% Specify which of the runs (1, ..., 10) to plot.
% Figures in JIN use 1st run, change this to see other runs
plot_run=1;  % Value should be in range 1 - 10.

figure
for ind=1:n_set
  ch_max=ch_max_all(ind,1);
  subplot(subplot_1,subplot_2,ind)
  plot(f12(:,1,plot_run,ind),f12(:,4,plot_run,ind),'k','LineWidth',line_width)
  hold on
  plot(f12(:,1,plot_run,ind),f12(:,cl12(plot_run,ind).col_R2p,plot_run,ind),'r','LineWidth',line_width)
  plot(f12(:,1,plot_run,ind),f12(:,cl12(plot_run,ind).col_R2n,plot_run,ind),'c','LineWidth',line_width)
  plot([0 freq_all(ind)],[cl12(plot_run,ind).ch_c95 cl12(plot_run,ind).ch_c95],'k--')
  axis([0 freq_all(ind) 0 ch_max+1e-7])
  box off
  hold off
  if (ind>freq_label_ind)
    xlabel('Frequency (Hz)')
  end
  set(gca,'TickDir','out')
  text(0.08*freq_all(ind),ch_max,[label_text(ind),')'],'FontSize',10)
end

%-----------------------------------------------------------------------------
% Plotting for neurons 1->3, frequency domain, Figure 3 in JIN
freq_all=[250 250 250 250 250 100 250 250];

figure
for ind=1:n_set
  ch_max=ch_max_all(ind,2);
  subplot(subplot_1,subplot_2,ind)
  plot(f13(:,1,plot_run,ind),f13(:,4,plot_run,ind),'k','LineWidth',line_width)
  hold on
  plot(f13(:,1,plot_run,ind),f13(:,cl13(plot_run,ind).col_R2p,plot_run,ind),'r','LineWidth',line_width)
  plot(f13(:,1,plot_run,ind),f13(:,cl13(plot_run,ind).col_R2n,plot_run,ind),'c','LineWidth',line_width)
  plot([0 freq_all(ind)],[cl13(plot_run,ind).ch_c95 cl13(plot_run,ind).ch_c95],'k--')
  axis([0 freq_all(ind) 0 ch_max])
  box off
  hold off
  if (ind>freq_label_ind)
    xlabel('Frequency (Hz)')
  end 
  set(gca,'TickDir','out')
  text(0.08*freq_all(ind),ch_max,[label_text(ind),')'],'FontSize',10)
end

%-----------------------------------------------------------------------------
% Plotting for neurons 1->3, time domain, Figure 4 in JIN
lag_range_all_13=[-10 40;-40 10;-50 50;-50 50;-30 30;-50 50;-10 75;-75 30];
rho_axis=[-0.02 0.14];
rho_axis_6=[-0.025 0.025];
label_text=['a','b','c','d','e','f','g','h'];

% Landscape format
lag_label_ind=4;
line_width=0.5;
subplot_1=2;
subplot_2=4;

figure
for ind=1:n_set
  lag_range=lag_range_all_13(ind,:);
  subplot(subplot_1,subplot_2,ind)
  plot(t13(:,1,plot_run,ind),t13(:,3,plot_run,ind),'k','LineWidth',line_width)
  hold on
  plot(lag_range,[0 0],'k--')
  plot(lag_range,+cl13(plot_run,ind).rho_c95*[1 1],'k-')
  plot(lag_range,-cl13(plot_run,ind).rho_c95*[1 1],'k-')
  plot([0 0],rho_axis,'k:')
  if ind==6
    axis([lag_range rho_axis_6])
    text(lag_range(1)+0.08*(lag_range(2)-lag_range(1)),rho_axis_6(2),[label_text(ind),')'],'FontSize',10)
  else
    axis([lag_range rho_axis])
    text(lag_range(1)+0.08*(lag_range(2)-lag_range(1)),rho_axis(2),[label_text(ind),')'],'FontSize',10)
  end 
  box off
  hold off
  if (ind>lag_label_ind)
    xlabel('Lag (ms)')
  end 
end

%-----------------------------------------------------------------------------
% Additional plots for individual configurations using psp2_R2 plotting
%  psp2_R2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)

% List of configurations (range 1 - 8)
plot_config_list=[1 2 3 4 5 8 8];

% List of runs (range 1 - 10) - should be same length as above
plot_run_list=[1 1 1 1 1 1 10];
% Edit these arrays to see other combinations, should have same number of entries in each.

%Plotting parameters
freq=200;
lag_tot=150;
lag_neg=75;
ch_max=0.6;
   
for ind=1:length(plot_config_list)
  run_no=plot_run_list(ind);
  config_no=plot_config_list(ind);
  figure
  psp2_R2(f12(:,:,run_no,config_no),t12(:,:,run_no,config_no),cl12(run_no,config_no),freq,lag_tot,lag_neg,ch_max)
  figure
  psp2_R2(f13(:,:,run_no,config_no),t13(:,:,run_no,config_no),cl13(run_no,config_no),freq,lag_tot,lag_neg,ch_max)
end
   