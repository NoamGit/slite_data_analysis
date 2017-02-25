% Script: R2_cn_pc_demo1.m
% Script to demonstrate use of NeuroSpec 2.11
% R2 conditional directionality analysis on spike train data from simulated Cortical Neuron networks.
%
% Copyright 2016, David M. Halliday.
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

% This script uses simulated cortical neuron spike train data, generates Figures 2 and 3 in JNM article:
%  Halliday DM, Senik MH, Stevenson CW, Mason R (2016). Non-parametric directionality analysis – Extension
%  for removal of a single common predictor and application to time series.
%  Journal of Neuroscience Methods, 268, 87–97. doi:10.1016/j.jneumeth.2016.05.008

% Load data
load JNM_dat1
% Contents of MAT file:
%  cn01_3pc_a, cn01_3pc_b, cn01_3pc_c, cn01_3pc_d:  Data matrices for configuration a), b), c), d), respectively.
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
% Storage for results from neurons 2-3: unconditional (coherence) and conditional (partial coherence)
n_set=length(filstem_list);
% R2 frequency domain matrices - 12 columns
f23=zeros(seg_size/2,12,n_runs,n_set);
f23p=zeros(seg_size/2,12,n_runs,n_set);
% R2 time domain matrices - 3 columns
t23=zeros(seg_size,3,n_runs,n_set);
t23p=zeros(seg_size,3,n_runs,n_set);

% R2 analysis of neurons 2-3
for set_no=1:n_set
  fil_stem=filstem_list{set_no};
  % select single set of data to use
  dat=eval(fil_stem);
  disp(' ')
  disp(['File no: ',num2str(set_no),' - ',fil_stem])

  % R2 analysis for individual runs
  for ind_run=1:n_runs
    % Neurons 2-3, unconditional analysis
    [f23(:,:,ind_run,set_no),t23(:,:,ind_run,set_no),cl23(ind_run,set_no)]=sp2a2_R2(dat(:,2,ind_run),dat(:,3,ind_run),rate,seg_pwr);
    cl23(ind_run,set_no).what=[fil_stem,': 2-3, run: ',num2str(ind_run)];
    cl23(ind_run,set_no).N1=sp_tot_all(ind_run,2,set_no);  % Append spike counts
    cl23(ind_run,set_no).N2=sp_tot_all(ind_run,3,set_no);
    cl23(ind_run,set_no).P1=sp_Pn_all(ind_run,2,set_no);   % Append expected values for Point Process spectra
    cl23(ind_run,set_no).P2=sp_Pn_all(ind_run,3,set_no);
    
    % Neurons 2-3|1, conditional, partial coherence analysis
    [f23p(:,:,ind_run,set_no),t23p(:,:,ind_run,set_no),cl23p(ind_run,set_no)]=sp2a2_R2_pc1(dat(:,2,ind_run),dat(:,3,ind_run),dat(:,1,ind_run),rate,seg_pwr);
    cl23p(ind_run,set_no).what=[fil_stem,': 2-3/1, run: ',num2str(ind_run)];
    cl23p(ind_run,set_no).N1=sp_tot_all(ind_run,2,set_no);  % Append spike counts
    cl23p(ind_run,set_no).N2=sp_tot_all(ind_run,3,set_no);
    cl23p(ind_run,set_no).P1=sp_Pn_all(ind_run,2,set_no);   % Append expected values for Point Process spectra
    cl23p(ind_run,set_no).P2=sp_Pn_all(ind_run,3,set_no);
  end 
end

%-----------------------------------------------------------------------------
% Plotting for neurons 2->3, unconditional analysis in time domain, Figure 2 in JNM
lag_range_all_23=[-50 50;-50 50;-50 50;-50 50];
rho_axis=[-0.021 0.035];
label_text=['a','b','c','d'];

% Landscape format
lag_label_ind=2;
line_width=0.5;
subplot_1=2;
subplot_2=2;

% Specify which of the runs (1, ..., 10) to plot.
% Figures in JIN use 1st run, change this to see other runs
plot_run=1;  % Value should be in range 1 - 10.

figure
for ind=1:n_set
  lag_range=lag_range_all_23(ind,:);
  subplot(subplot_1,subplot_2,ind)
  plot(t23(:,1,plot_run,ind),t23(:,3,plot_run,ind),'k','LineWidth',line_width)
  hold on
  plot(lag_range,[0 0],'k--')
  plot(lag_range,+cl23(plot_run,ind).rho_c95*[1 1],'k-')
  plot(lag_range,-cl23(plot_run,ind).rho_c95*[1 1],'k-')
  plot([0 0],rho_axis,'k:')
  axis([lag_range rho_axis])
  text(lag_range(1)+0.08*(lag_range(2)-lag_range(1)),rho_axis(2),[label_text(ind),')'],'FontSize',10)
  box off
  hold off
  if (ind>lag_label_ind)
    xlabel('Lag (ms)')
  end 
end

%-----------------------------------------------------------------------------
% Plotting for neurons 2->3/1, conditional analysis in time domain, Figure 3 in JNM
% Using same plot parameters as above

figure
for ind=1:n_set
  lag_range=lag_range_all_23(ind,:);
  subplot(subplot_1,subplot_2,ind)
  plot(t23p(:,1,plot_run,ind),t23p(:,3,plot_run,ind),'k','LineWidth',line_width)
  hold on
  plot(lag_range,[0 0],'k--')
  plot(lag_range,+cl23p(plot_run,ind).rho_c95*[1 1],'k-')
  plot(lag_range,-cl23p(plot_run,ind).rho_c95*[1 1],'k-')
  plot([0 0],rho_axis,'k:')
  axis([lag_range rho_axis])
  text(lag_range(1)+0.08*(lag_range(2)-lag_range(1)),rho_axis(2),[label_text(ind),')'],'FontSize',10)
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
plot_config_list=[1 2 3 4];

% List of runs (range 1 - 10) - should be same length as above
% plot_run_list=[1 1 1 1];  % first run of each set
plot_run_list=[3 5 6 9];    % arbitrary selection
% Edit these arrays to see other combinations, should have same number of entries in each.

%Plotting parameters
freq=200;
lag_tot=150;
lag_neg=75;
ch_max=0.15;
   
for ind=1:length(plot_config_list)
  run_no=plot_run_list(ind);
  config_no=plot_config_list(ind);
  figure
  psp2_R2(f23(:,:,run_no,config_no),t23(:,:,run_no,config_no),cl23(run_no,config_no),freq,lag_tot,lag_neg,ch_max)
  figure
  psp2_R2(f23p(:,:,run_no,config_no),t23p(:,:,run_no,config_no),cl23p(run_no,config_no),freq,lag_tot,lag_neg,ch_max)
end
   