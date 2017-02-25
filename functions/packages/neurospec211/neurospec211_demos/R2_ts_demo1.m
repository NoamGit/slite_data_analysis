% Script: R2_ts_demo1.m
% Script to demonstrate use of NeuroSpec 2.1
% R2 directionality analysis on coupled bivariate time series data
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

% Generate uncorrelated bivariate data
% Either white noise as N(0,1) or coloured noise as AR(1).
%
% Select noise_type=0 OR noise_type=1
noise_type=0;   % Use this for white noise, N(0,1) processes
%noise_type=1;  % Use this for coloured noise, AR(1) processes

pts_tot=100000;  % No of data points to generate, 10^5
% Storage for data matrix - 5 columns.
dat=zeros(pts_tot,5);

switch noise_type
  case 0
    % Normally distributed, zero mean, unit variance: N(0,1).
    dat=randn(pts_tot,2);
  case 1  
    % Non white signals, generated as AR(1) processes using external function AR1_gen1
    dat(:,1)=AR1_gen1(pts_tot);
    dat(:,2)=AR1_gen1(pts_tot);
  otherwise
    return
end

% Add third column which is average of first two.
dat(:,3)=0.5*(dat(:,1)+dat(:,2));

% Column 4 is delayed version of column 3. Delay of n_del samples.
n_del=30;  % Change this for different delay
dat(1:n_del,4)=randn(n_del,1);                 % Initial values randomly generated
dat(n_del+1:pts_tot,4)=dat(1:pts_tot-n_del,3); % Rest copied from column 3

% Generate filtered version of column 1.
% Using 4 term Kaiser-Bessel window, length n_filt
n_filt=15;
% Fixed KB coefficients, from Harris 1978, dx.doi.org/10.1109/PROC.1978.10837
a0=0.40243;
a1=0.49804;
a2=0.09831;
a3=0.00122;
% Window is w(n) = a0 - a1 cos(2 pi n/N) + a2 cos (2 pi 2n/N) - a3 cos(2 pi 3n/N), n=0:N-1 (Harris, 1978, P65)
a=1;
ind_w=2*pi*(0:n_filt-1)'/n_filt;
b=a0-a1*cos(ind_w)+a2*cos(2*ind_w)-a3*cos(3*ind_w);
b=b/sum(b);
dat_filt=filter(b,a,dat(:,1));
% Column 5 is average of filtered version of column 1 and column 2.
dat(:,5)=0.5*(dat_filt+dat(:,2));
% Now have 5 column data matrix to use for directionality analysis

% Not needed for analysis
clear a a0 a1 a2 a3 b ind_w dat_filt

%------------------------------------------------------------------------
% For purposes of analysis: Assume sampling rate is 1000/sec
samp_rate=1000;

% Set power for segment length as 10.
% Thus T=2^10 = 1024.
% 1024 points @ 1000/sec sampling gives segment length of 1.024 sec
% Frequency resolution is inverse of this 0.977 Hz.
seg_pwr=10;

%------------------------------------------------------------------------
% R2 analysis on columns 1 - 3, instantaneous (zero lag) correlation
%  Analysis:  [f,t,cl]=sp2a2_R2(x,y,samp_rate,seg_pwr) 
%  Plotting:  psp2_R2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)

[f1,t1,cl1]=sp2a2_R2(dat(:,1),dat(:,3),samp_rate,seg_pwr);
cl1.what='column 1 - 3';
figure
% Plotting, freq: 500; lags:100,50; ch_max=1
psp2_R2(f1,t1,cl1,500,100,50,1)

%------------------------------------------------------------------------
% R2 analysis on columns 2 - 4, correlation in single time bin with delay n_del

[f2,t2,cl2]=sp2a2_R2(dat(:,2),dat(:,4),samp_rate,seg_pwr);
cl2.what=['Column 2 - 4. Delay: ',num2str(n_del)];
figure
% Plotting, freq: 500; lags:100,50; ch_max=1
psp2_R2(f2,t2,cl2,500,100,50,1)

%------------------------------------------------------------------------
% R2 analysis on columns 1 - 5, correlation through KB filter of length n_filt

[f3,t3,cl3]=sp2a2_R2(dat(:,1),dat(:,5),samp_rate,seg_pwr);
cl3.what=['Column 1 - 5. Filter: ',num2str(n_filt),' pts.'];
figure
% Plotting, freq: 200; lags:100,50; ch_max=1
psp2_R2(f3,t3,cl3,200,100,50,1)
