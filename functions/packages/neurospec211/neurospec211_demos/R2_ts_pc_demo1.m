% Script: R2_ts_pc_demo1.m
% Script to demonstrate use of NeuroSpec 2.11
% R2 conditional directionality analysis on delayed mixtures of Gausssian signals
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

% Conditional directionality analysis on delayed mixtures of Gausssian signals, generates Figure 4 in JNM article:
%  Halliday DM, Senik MH, Stevenson CW, Mason R (2016). Non-parametric directionality analysis – Extension
%  for removal of a single common predictor and application to time series.
%  Journal of Neuroscience Methods, 268, 87–97. doi:10.1016/j.jneumeth.2016.05.008
%
% Model generated data: incorporates delays to induce directionality between x(t) and y(t)
%  x(t) =  a1*z1(t-1) + a2*z2(t)   + sqrt(1-(a1*a1+a2*a2))*e1(t)
%  y(t) =  a1*z1(t)   + a2*z2(t-1) + sqrt(1-(a1*a1+a2*a2))*e2(t)
% All inputs have unit variance: var{z1}=var(z2)=var{e1}=var{e2} = 1
%
% Parameters are weigths: a1, a2
%
% Ordinary coherence: |R^2_yx(w)|^2  =  a1^4 + a2^4 + 2 * a1^2 * a2^2 * cos(2 w tau)
%  w is radian frequency, tau is delay (default below is one time step)
%
% Partial coherence  R^2_xy/z1  =  (a2)^2 / (1-a1^2)^2. This will be in forward direction x -> y
%                    R^2_xy/z2  =  (a1)^2 / (1-a2^2)^2. This will be in reverse direction y -> x
%
% Forward and reverse components for unconditional directional analysis, decomposition of ordinary coherence
%  Forward direction, ratio for |R'_yx;+|^2 = a2^4 / (a1^4 + a2^4)
%  Reverse direction, ratio for |R'_yx;-|^2 = a1^4 / (a1^4 + a2^4)

% Specify target coherence and fractional contribution from a1
r2=0.9;        % Target coherence at zero frequency
a1_frac=2/3;   % Fractional contribution from z1(t) through weight a1

% Number of data points to generate
pts_tot=10^5;
% Specify delay as number of time steps. Default is 1 time step
del_pts=1; % No of time steps

% Calculate weights a1, a2
a1=sqrt(sqrt(r2)*a1_frac);
a2=sqrt(sqrt(r2)*(1-a1_frac));

% Generate data matrices with 4 random processes: z1, z2; e1, e2
dat_z=randn(pts_tot+del_pts,2);
dat_e=randn(pts_tot+del_pts,2);
% Shouldn't need but subtract mean and standardize variance
dat_z(:,1)=(dat_z(:,1)-mean(dat_z(:,1)))/std(dat_z(:,1));  % z1(t)
dat_z(:,2)=(dat_z(:,2)-mean(dat_z(:,2)))/std(dat_z(:,2));  % z2(t)
dat_e(:,1)=(dat_e(:,1)-mean(dat_e(:,1)))/std(dat_e(:,1));  % e1(t)
dat_e(:,2)=(dat_e(:,2)-mean(dat_e(:,2)))/std(dat_e(:,2));  % e2(t)

% Generate x(t)
dat_x=a1*dat_z(1:pts_tot,1)+a2*dat_z(1+del_pts:pts_tot+del_pts,2)+sqrt(1-(a1*a1+a2*a2))*dat_e(1:pts_tot,1);
% Generate y(t)
dat_y=a1*dat_z(1+del_pts:pts_tot+del_pts,1)+a2*dat_z(1:pts_tot,2)+sqrt(1-(a1*a1+a2*a2))*dat_e(1:pts_tot,2);

% Combined data matrix for analysis
dat=[dat_z(1:pts_tot,:) dat_x dat_y];

%--------------------------------------
% Directionality analysis
% Use an assumed sampling rate of 1000/sec
rate=1000;
delta_t=1/rate;  % dt
% segment size for analysis T=2^10
seg_pwr=10;

% Unconditional directionality analysis
[f34,t34,cl34]=sp2a2_R2(dat(:,3),dat(:,4),rate,seg_pwr);
cl34.what=['Ch: 3-4. Target R2: ',num2str([(a1*a1+a2*a2)^2]),'. Del: ',num2str(del_pts)];

% Conditional directionality analysis with z1(t) and z2(t) as predictors
[f341p,t341p,cl341p]=sp2a2_R2_pc1(dat(:,3),dat(:,4),dat(:,1),rate,seg_pwr);  % z1(t) as predictor
[f342p,t342p,cl342p]=sp2a2_R2_pc1(dat(:,3),dat(:,4),dat(:,2),rate,seg_pwr);  % z2(t) as predictor
cl341p.what=['Partial: 3-4/1. Target R2: ',num2str(a2^4/(1-a1*a1)^2),'. Del: ',num2str(del_pts)];
cl342p.what=['Partial: 3-4/2. Target R2: ',num2str(a1^4/(1-a2*a2)^2),'. Del: ',num2str(del_pts)];

% Theoretical ordinary coherence
r2_del=(a1*a1)^2 + (a2*a2)^2 + 2*a1*a1*a2*a2*cos(2*f34(:,1)*2*pi*delta_t*del_pts);
% Directional components for ordinary coherence
r2_ratio_f=a2^4/(a1^4+a2^4); % ratio for forward component, |R'_yx;+|^2
r2_ratio_r=a1^4/(a1^4+a2^4); % ratio for reverse component, |R'_yx;-|^2

% Calculate theoretical overall R^2_yx value by integration of theoretical coherence
R2yx=mean(r2_del);

% Conditional analysis - theoretical partial coherence values
r2_pc_z1=a2^4/(1-a1*a1)^2;  % With z1 as predictor. Constant partial coherence, in forward direction
r2_pc_z2=a1^4/(1-a2*a2)^2;  % With z2 as predictor. Constant partial coherence, in reverse direction

%--------------------------------------
% Plotting -  unconditional analysis
figure
psp2_R2(f34,t34,cl34,500,20,10,1)
% Add in theoretical coherence
subplot(4,2,3)
hold on
plot(f34(:,1),r2_del,'r','LineWidth',2)
hold off

% Include enlarged directional components with target values superimposed
subplot(2,2,3)
psp_ch1_R2(f34,cl34,500,1)
hold on
plot(f34(:,1),r2_ratio_f*r2_del,'b','LineWidth',2) % Forward direction target in blue on red estimate
plot(f34(:,1),r2_ratio_r*r2_del,'r','LineWidth',2) % Reverse direction target in red on blue estimate

% Display theoretical values
disp(['Theoretical - ordinary R^2: ',num2str([R2yx]),',   Partial R^2: ',num2str([a2^4/(1-a1*a1)^2  a1^4/(1-a2*a2)^2])])

%--------------------------------------
% Plotting -  conditional analysis, partial coherence with z1 as predictor
figure
psp2_R2(f341p,t341p,cl341p,500,20,10,1)
% Add in target partial coherence, |R_yx/z1|^2. This forward, plot in red
subplot(4,2,3)
hold on
plot(f34(:,1),r2_pc_z1*ones(length(f341p),1),'r','LineWidth',2)
hold off

%--------------------------------------
% Plotting -  conditional analysis, partial coherence with z2 as predictor
figure
psp2_R2(f342p,t342p,cl342p,500,20,10,1)
% Add in target partial coherence, |R_yx/z2|^2. This forward, plot in blue
subplot(4,2,3)
hold on
plot(f34(:,1),a1^4/(1-a2*a2)^2*ones(length(f342p),1),'b','LineWidth',2)

disp(' ')

%----------------------------------------------------------
% JNM  - Figure 4.
% Multiply frequency values by 1/rate to get normalised frequency scale, where fN = 0.5
f_fac=1/rate;
figure
subplot(2,2,1)
plot(f_fac*f34(:,1),f34(:,4),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on
plot(f_fac*f34(:,1),r2_del,'k','LineWidth',1)
box off
hold off
axis([0 0.5 0 1])
text(0.03,1,'a)','FontSize',10)
text('Interpreter','latex','String','$$|\hat{R}_{yx}(\lambda)|^2$$','Position',[.15 .7],'FontSize',10)

subplot(2,2,3)
plot(f_fac*f34(:,1),f34(:,cl34.col_R2p),'Color',[0 0 0],'LineWidth',1)
hold on
plot(f_fac*f34(:,1),f34(:,cl34.col_R2n),'Color',[0.5 0.5 0.5],'LineWidth',1)
% Now adding theoretical values
plot(f_fac*f34(:,1),r2_ratio_r*r2_del,'Color',[0 0 0],'LineWidth',1)
plot(f_fac*f34(:,1),r2_ratio_f*r2_del,'Color',[0.5 0.5 0.5],'LineWidth',1)
box off
hold off
xlabel('Frequency')
axis([0 0.5 0 1])
text(0.03,1,'c)','FontSize',10)
text('Interpreter','latex','String','$$|\hat{R}''_{yx;-}(\lambda)|^2$$','Position',[.1,.97],'FontSize',10)
text('Interpreter','latex','String','$$|\hat{R}''_{yx;+}(\lambda)|^2$$','Position',[.1 .77],'FontSize',10)
% Lines
hold on
plot([0.4 0.45],[0.97 0.97],'Color',[0.5 0.5 0.5],'LineWidth',1)
plot([0.4 0.45],[0.77 0.77],'Color',[0   0   0  ],'LineWidth',1)


subplot(2,2,2)
plot(f_fac*f341p(:,1),f341p(:,cl341p.col_R2p),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on
plot(f_fac*f341p(:,1),r2_pc_z1*ones(length(f341p),1),'k','LineWidth',1)
box off
hold off
axis([0 0.5 0 1])
text(0.03,1,'b)','FontSize',10)
text('Interpreter','latex','String','$$|\hat{R}''_{yx|z_1;+}(\lambda)|^2$$','Position',[.1 .45],'FontSize',10)

subplot(2,2,4)
plot(f_fac*f342p(:,1),f342p(:,cl342p.col_R2n),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on
plot(f_fac*f342p(:,1),r2_pc_z2*ones(length(f342p),1),'k','LineWidth',1)
box off
hold off
axis([0 0.5 0 1])
xlabel('Frequency')
text(0.03,1,'d)','FontSize',10)
text('Interpreter','latex','String','$$|\hat{R}''_{yx|z_2;-}(\lambda)|^2$$','Position',[.1 .55],'FontSize',10)

