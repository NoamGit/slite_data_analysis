% Script: R2_pc_table1_demo.m
% Script to demonstrate use of NeuroSpec 2.11
% R2 conditional directionality analysis on mixtures of Gausssian signals
% Generates Table 1 in JNM paper (reference below)
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

% Conditional directionality analysis on mixtures of Gausssian signals, generates Table 1 in JNM paper:
%  Halliday DM, Senik MH, Stevenson CW, Mason R (2016). Non-parametric directionality analysis – Extension
%  for removal of a single common predictor and application to time series.
%  Journal of Neuroscience Methods, 268, 87–97. doi:10.1016/j.jneumeth.2016.05.008
%
% Model generated data: two random processes with single common input
%  x(t) =  a1*z(t) + sqrt(1-(a1*a1))*e1(t)
%  y(t) =  a1*z(t) + sqrt(1-(a1*a1))*e2(t)
%     where            (a1*a1) < 1  
%
% All inputs have unit variance: var{z}=var{e1}=var{e2} = 1
%
% Parameter is weight: a1

% Ordinary coherence R^2_yx     =  (a1^2)^2
% Partial coherence  R^2_yx/z   =  zero

% Specify target R^2 between x and y
target_r12_all=[0.1 0.3 0.5 0.7 0.9];

% No of runs at each Target R^2
run_tot=100;
% No of segments for each run
seg_tot=100;

% Analysis parameters
rate=1000;                 % Assumed sampling rate of 1000/sec
seg_pwr=10;                % Segment power
seg_size=2^seg_pwr;        % T=2^10
pts_tot=seg_size*seg_tot;  % Number of data points

% Storage for metrics
R2_yx_all=zeros(run_tot,length(target_r12_all));
R2_yxz_all=zeros(run_tot,length(target_r12_all));

for ind_r2=1:length(target_r12_all)

  target_r12=target_r12_all(ind_r2);
  disp(' ')
  disp(['Target R^2: ',num2str(target_r12)])
  disp(' ')

  for ind_run=1:run_tot

    % Generates 3 random processes z, e1, e2
    dat_n=randn(pts_tot,3);
    % Shouldn't need but standardize anyway
    dat_n(:,1)=(dat_n(:,1)-mean(dat_n(:,1)))/std(dat_n(:,1));  % z1(t)
    dat_n(:,2)=(dat_n(:,2)-mean(dat_n(:,2)))/std(dat_n(:,2));  % e1(t)
    dat_n(:,3)=(dat_n(:,3)-mean(dat_n(:,3)))/std(dat_n(:,3));  % e2(t)
    
    % Generate x(t) and y(t)    
    a1=sqrt(sqrt(target_r12));                       
    dat_x=a1*dat_n(:,1)+sqrt(1-a1*a1)*dat_n(:,2); % x(t)
    dat_y=a1*dat_n(:,1)+sqrt(1-a1*a1)*dat_n(:,3); % y(t)
    
    % Unconditional analysis between x and y: R^2_yx
    [f23,t23,cl23]=sp2a2_R2(dat_x,dat_y,rate,seg_pwr);

    % Conditional analysis between x and y with z as predictor: R^2_yx|z
    [f23p,t23p,cl23p]=sp2a2_R2_pc1(dat_x,dat_y,dat_n(:,1),rate,seg_pwr);

    % Store R^2 for each run
    R2_yx_all(ind_run,ind_r2)=cl23.R2_ch;
    R2_yxz_all(ind_run,ind_r2)=cl23p.R2_ch;
  
  end
  
end

% Histograms of R^2 values for unconditional analysis
figure
hist_bin=20;
for plot_ind=1:length(target_r12_all)
  subplot(length(target_r12_all),1,plot_ind)
  hist(R2_yx_all(:,plot_ind),hist_bin)
  title(['Target R^2: ',num2str(target_r12_all(plot_ind))])
end

% Results
disp(' ')
disp(' ')
disp(['R^2_yx    Mean: ',num2str(mean(R2_yx_all),3)])
disp(['R^2_yx     Std: ',num2str(std(R2_yx_all),3)])
disp(['R^2_yx|z  Mean: ',num2str(mean(R2_yxz_all),3)])
disp(['R^2_yx|z   Std: ',num2str(std(R2_yxz_all),3)])

% Generate Table 1 in paper
Table_1=[
target_r12_all
mean(R2_yx_all)
mean(R2_yx_all)-2*std(R2_yx_all)
mean(R2_yx_all)+2*std(R2_yx_all)
mean(R2_yxz_all)
mean(R2_yxz_all)-2*std(R2_yxz_all)
mean(R2_yxz_all)+2*std(R2_yxz_all)
]';
 
% Display table
Table_1
