% Script: R2_AR4_demo1.m
% Script to demonstrate use of NeuroSpec 2.1
% R2 directionality analysis on coupled bivariate AR(4) data
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

% Generate correlated bivariate AR(4) data using external function AR4_gen1, with 10^5 samples
dat=AR4_gen1(100000);

% For analysis - set nominal sampling rate of 1, Nyquist fN=0.5
rate=1;
seg_pwr=10;
% R2 analysis using average periodogram
[f1,t1,cl1]=sp2a2_R2(dat(:,1),dat(:,2),rate,seg_pwr);
cl1.what='AR(4), correlated, 1->2';

% R2 analysis using multi-taper, NW=2.
[f2,t2,cl2]=sp2a2_R2_mt(dat(:,1),dat(:,2),rate,seg_pwr,'M2');
cl2.what='AR(4), correlated, 1->2 (MT)';

%---------------------------------------------
% Plotting parameters
freq=0.2;         % Fractional frequency, where fN=0.5
lag_tot=300*1000; % Lags specified in ms, his is 300 (one second) bins in ms
lag_neg=100*1000; % Negative lags over 100 bins
ch_max=1;         % Max coherence

% Plotting spectral and directionality analysis
figure
psp2_R2(f1,t1,cl1,freq,lag_tot,lag_neg,ch_max)
figure
psp2_R2(f2,t2,cl2,freq,lag_tot,lag_neg,ch_max)

% R2 metrics over restricted frequency range, use same range as plotting
cl1_f=R2_f1(f1,cl1,freq);
cl2_f=R2_f1(f2,cl2,freq);

% Display metrics and percentage of R2 in forward direction
disp(' ')
disp(['R2 (0 p n): ',num2str(cl1_f.R2_ch_f),' (',num2str([cl1_f.R2_0_f cl1_f.R2_p_f cl1_f.R2_n_f]),')'])
disp(['R2 percentage forward: ',num2str(cl1_f.R2_per_p)])
disp(' ')
disp(['MT R2 (0 p n): ',num2str(cl2_f.R2_ch_f),' (',num2str([cl2_f.R2_0_f cl2_f.R2_p_f cl2_f.R2_n_f]),')'])
disp(['MT R2 percentage forward: ',num2str(cl2_f.R2_per_p)])

