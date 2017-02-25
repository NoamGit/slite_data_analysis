function [f,t,cl]=sp2a2_R2_pc1(x,y,z,samp_rate,seg_pwr)
% function [f,t,cl]=sp2a2_R2_pc1(x,y,z,samp_rate,seg_pwr)
%
% Three channel partial spectra analysis, with conditional non-parametric directionality analysis.
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
%
% Input arguments
%  x          Channel 1 (input)     time series vector
%  y          Channel 2 (output)    time series vector
%  z          Channel 3 (predictor) time series vector
%  samp_rate  Sampling rate (samples/sec)
%  seg_pwr    Segment length - specified as power of 2
%
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates
%  t    column matrix with      time domain parameter estimates
%  cl   single structure with scalar values related to analysis
%
% Output parameters
%  f column 1       frequency in Hz
%  f column 2       Log input/x  partial spectrum
%  f column 3       Log output/y partial spectrum
%  f column 4       Partial coherence
%  f column 5       Partial phase
%  f column 6       MMSE whitened partial spectrum for process x, JNM (16)
%  f column 7       MMSE whitened partial spectrum for process y, JNM (16)
%  f column 8       Partial coherence between whitened processes, JNM (17)
%  f column 9       Partial phase between whitened processes
%  f column 10      Zero lag partial coherence component, JNM (28)
%  f column 11      Forward  partial coherence component, JNM (29)
%  f column 12      Reverse  partial coherence component, JNM (27)
%
%  t column 1       Lag in ms.
%  t column 2       Partial cumulant/covariance density.
%  t column 3       rhoyxz, JNM (19)
%
%  cl.type          Analysis type 0, for compatability with NeuroSpec2.0
%  cl.seg_size      Segment length, T
%  cl.seg_tot       Number of segments, L
%  cl.seg_tot_var   Effective no of segments, L', used to calculate confidence limits
%  cl.samp_tot      Number of samples analysed, R=LT
%  cl.samp_rate     Sampling rate of data (samples/sec)
%  cl.dt            Time domain bin width (ms)
%  cl.df            Frequency domain bin width (Hz)
%  cl.f_c95         95% confidence limit for Log partial spectral estimates
%  cl.ch_c95        95% confidence limit for partial coherence estimate
%  cl.q_c95         95% confidence limit for partial cumulant density estimate
%  cl.rho_c95       95% confidence limit for rho estimate
%  cl.col_R20       column in f matrix containing zero lag partial coherence component
%  cl.col_R2p       column in f matrix containing forward (positive lag) partial coherence component
%  cl.col_R2n       column in f matrix containing reverse (negative lag) partial coherence component
%  cl.col_rho       column in t matrix containing rho estimate
%  cl.R2            R2 value, R^2_yx|z, JNM (22)
%  cl.R2_0          Component of R2 at zero lag, R^2_yx|z;0
%  cl.R2_p          Component of R2 in forward direction, positive lag component, R^2_yx|z;+
%  cl.R2_n          Component of R2 in reverse direction, negative lag component, R^2_yx|z;-
%  cl.R2_ch         R2 estimated from partial coherence between original processes
%  cl.R2_chw        R2 estimated from partial coherence between MMSE whitened processes
%  cl.R2_fprime_0   Component of R2 at zero lag estimated from f'
%  cl.R2_fprime_p   Component of R2 at positive lag estimated from f'
%  cl.R2_fprime_n   Component of R2 at negative lag estimated from f'
%  cl.N1, cl.N2     Spike counts for point process data
%  cl.P1, cl.P2     Spike rates for point process data
%  cl.opt_str       Copy of options string, empty here
%  cl.what          Text label, used in plotting routines
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%     Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Halliday DM (2015) Nonparametric directionality measures for time series and point process data,
%     Journal of Integrative Neuroscience, 14(2), 253-277 DOI: 10.1142/S0219635215300127
% 3. Halliday DM, Senik MH, Stevenson CW, Mason R (2016) Non-parametric directionality analysis – Extension
%     for removal of a single common predictor and application to time series.
%     Journal of Neuroscience Methods 268, 87-97. DOI: 10.1016/j.jneumeth.2016.05.008
%
% These references referred to in comments as PBMB, JIN, JNM respectively
%
% function [f,t,cl]=sp2a2_R2_pc1(x,y,z,samp_rate,seg_pwr)

% Check number of input arguments
if (nargin<5)
  error(' Not enough input arguments');
end
% Check number of output arguments
if (nargout<3)
  error(' Not enough output arguments');
end

% Check for single column data
[nrow,ncol]=size(x);
if (ncol~=1)
  error(' Input NOT single column: x')
end
[nrow,ncol]=size(y);
if (ncol~=1)
  error(' Input NOT single column: y')
end
[nrow,ncol]=size(z);
if (ncol~=1)
  error(' Input NOT single column: z')
end

pts_tot=length(x);           % Determine size of data vector
if (length(y)~=pts_tot)      % Check that input vectors are equal length
  error (' Unequal length data arrays: x, y');
end
if (length(z)~=pts_tot)      % Check that input vectors are equal length
  error (' Unequal length data arrays: x, z');
end

if (max(size(samp_rate)) ~= 1)
  error(' Non scalar value for: samp_rate');
end
if (max(size(seg_pwr)) ~= 1)
  error(' Non scalar value for: seg_pwr');
end

% Segment for FFT - Rows: seg_size, T. Columns: seg_tot, L
seg_size=2^seg_pwr;              % Segment length, T
seg_tot=fix(length(x)/seg_size); % No of segments, L
samp_tot=seg_tot*seg_size;       % Record length,  R=LT
x=reshape(x(1:samp_tot),seg_size,seg_tot);  % T rows, L columns
y=reshape(y(1:samp_tot),seg_size,seg_tot);  % T rows, L columns
z=reshape(z(1:samp_tot),seg_size,seg_tot);  % T rows, L columns

% Create zero mean sequence for each segment.
for ind=1:seg_tot
  x(:,ind)=x(:,ind)-mean(x(:,ind));
  y(:,ind)=y(:,ind)-mean(y(:,ind));
  z(:,ind)=z(:,ind)-mean(z(:,ind));
end

% FFT & average periodogram
dx=fft(x);   % dFT across columns, PBMB (4.1)
dy=fft(y);
dz=fft(z);

% Spectra fzz, fxz, fyz needed to generate conditional dFT
psd_fac=1/(2*pi*samp_tot);
fzz=psd_fac*sum(abs(dz.*dz),2);  % Using PBMB (5.2)
fxz=psd_fac*sum(dx.*conj(dz),2); % Using PBMB (5.2)
fyz=psd_fac*sum(dy.*conj(dz),2); % Using PBMB (5.2)

% Gain factors to generate conditional dFT, with zero frequency set to zero
gxz=[0; fxz(2:seg_size)./fzz(2:seg_size)];
gyz=[0; fyz(2:seg_size)./fzz(2:seg_size)];

% Generate conditional dFTs: dx|z, dy|z
dxz=dx;
dyz=dy;
for ind=1:seg_tot
  dxz(:,ind)=dxz(:,ind)-gxz.*dz(:,ind);  % JNM (8), JNM (30)
  dyz(:,ind)=dyz(:,ind)-gyz.*dz(:,ind);  % JNM (9), JNM (31)
end

% Estimate partial spectra from conditional dFT using average over segments
fxxz=psd_fac*sum(abs(dxz.*dxz),2);  % JNM (10), JNM (32)
fyyz=psd_fac*sum(abs(dyz.*dyz),2);  % JNM (11)
fyxz=psd_fac*sum(dyz.*conj(dxz),2); % Partial cross spectrum

% Estimate partial coherence JNM (3), PBMB (8.8), with zero frequency set to zero
chyxz=[0; abs(fyxz(2:seg_size).*fyxz(2:seg_size))./(fxxz(2:seg_size).*fyyz(2:seg_size))];

%Pre-whitening stage - generate MMSE filters, weight for zero frequency set to zero
wxxz=[0; 1./sqrt(fxxz(2:seg_size))];  % JNM (12), JNM (33)
wyyz=[0; 1./sqrt(fyyz(2:seg_size))];  % JNM (13), JNM (34)

% Apply MMSE pre-whitening filters to generate pre-whitened dFTs: dwx|z, dwy|z
dxzw=dxz;
dyzw=dyz;
for ind=1:seg_tot
  dxzw(:,ind)=dxzw(:,ind).*wxxz;  % JNM (14), JNM (35)
  dyzw(:,ind)=dyzw(:,ind).*wyyz;  % JNM (15), JNM (36)
end

% Re-calculate spectra using pre-whitened processes
fxxzw=psd_fac*sum(abs(dxzw.*dxzw),2);   % Spectra of conditioned MMSE whitened process x. Equals 1 at all freqs, JNM (16)
fyyzw=psd_fac*sum(abs(dyzw.*dyzw),2);   % Spectra of conditioned MMSE whitened process y. Equals 1 at all freqs, JNM (16)
fyxzw=psd_fac*sum(dyzw.*conj(dxzw),2);  % Complex cross spectrum between conditioned whitened processes
chyxzw=abs(fyxzw.*fyxzw);               % Partial coherence from fyxzw only, JNM (17), JNM (37)
rhoyxz=real(ifft(fyxzw));               % rhoyxz, no factors as 1/T in ifft(), estimate of JNM (19) 

% Estimate R2 in frequency domain
R2_ch=(1/seg_size)*sum(chyxz);   % Integral from -pi to +pi, JNM (5)
R2_chw=(1/seg_size)*sum(chyxzw); % Integral from -pi to +pi, JNM (18)

% Estimate R2 in time domain
R2_rho=sum(rhoyxz.^2);  % Integral over all lags, JNM (20)

% Estimate separate components, JNM (21), JNM (22)
R2_rho_0=rhoyxz(1).^2;                          % zero lag, bin 1,               R^2_yx|z;0
R2_rho_p=sum(rhoyxz(2:seg_size/2).^2);          % positive lags, bins 2:(T/2),   R^2_yx|z;+
R2_rho_n=sum(rhoyxz(seg_size/2+1:seg_size).^2); % negative lags, bins (T/2)+1:T, R^2_yx|z;-

% Decompose rho into three components by lag, to estimate f' functions
rhoyxz_3=zeros(length(rhoyxz),3);
rhoyxz_3(1,1)=rhoyxz(1);                                         % zero lag
rhoyxz_3(2:seg_size/2,2)=rhoyxz(2:seg_size/2);                   % positive lags
rhoyxz_3(seg_size/2+1:seg_size,3)=rhoyxz(seg_size/2+1:seg_size); % negative lags

% Switch back to frequency domain, these are f' estimates
f_prime=fft(rhoyxz_3);      % JNM (24, 25, 23): f'yx|z;0, f'yx|z;+, f'yx|z;-
f_prime2=abs(f_prime).^2;   % These are |f'|^2 values: |f'yx|z;0|^2, |f'yx|z;+|^2, |f'yx|z;-|^2

% Scale to get relative contributions to partial coherence from 0,+,- components
f_prime2_fac=sum(f_prime2,2);
R2_weight=f_prime2;
R2_weight(:,1)=R2_weight(:,1)./f_prime2_fac;  % Weighting factor in JNM (28), range [0 1]
R2_weight(:,2)=R2_weight(:,2)./f_prime2_fac;  % Weighting factor in JNM (29), range [0 1]
R2_weight(:,3)=R2_weight(:,3)./f_prime2_fac;  % Weighting factor in JNM (27), range [0 1]

% R^2 decomposition using integrated |f'|^2 - should agree with decomposition by lag
R2_fprime_0=(1/seg_size)*sum(f_prime2(:,1));  % These are terms in conditioned version of JIN (2.15)
R2_fprime_p=(1/seg_size)*sum(f_prime2(:,2));
R2_fprime_n=(1/seg_size)*sum(f_prime2(:,3));

%-----------------------------------------------------------------------
% Construct output matrix f
f=zeros(seg_size/2,12);
f_index=(2:seg_size/2+1)';    % Indexing for output, DC not output.
deltaf=samp_rate/seg_size;
f(:,1)=(f_index-1)*deltaf;    % Column 1 - frequencies in Hz
f(:,2)=log10(fxxz(f_index));  % Column 2 - Log partial spectrum ch 1, x, JNM (10)
f(:,3)=log10(fyyz(f_index));  % Column 3 - Log partial spectrum ch 2, y, JNM (11)
f(:,4)=chyxz(f_index);        % Column 4 - Partial coherence, PBMB (8.8), JNM (3)
f(:,5)=angle(fyxz(f_index));  % Column 5 - Partial phase, PBMB (8.9)
f(:,6)=fxxzw(f_index);        % Column 6 - MMSE whitened partial spectra process x, JNM (16)
f(:,7)=fyyzw(f_index);        % Column 7 - MMSE whitened partial spectra process y, JNM (16)
f(:,8)=chyxzw(f_index);       % Column 8 - Partial coherence between whitened processes, JNM (17)
f(:,9)=angle(fyxzw(f_index)); % Column 9 - Partial phase between whitened processes
f(:,10)=R2_weight(f_index,1).*chyxzw(f_index);   % Column 10 - Zero lag partial coherence component, JNM (28)
f(:,11)=R2_weight(f_index,2).*chyxzw(f_index);   % Column 11 - Forward partial  coherence component, JNM (29)
f(:,12)=R2_weight(f_index,3).*chyxzw(f_index);   % Column 12 - Reverse partial  coherence component, JNM (27)

% Construct time domain matrix t
% Time domain output: Col 1 - lag  range (-T/2)*dt to (T/2-1)*dt, Col 2 - partial cumulant, Col 3 - rho
deltat=1000.0/samp_rate;    % dt in msec.
t(:,1)=((1:seg_size)'-seg_size/2-1)*deltat;  % Time lag in ms
t([seg_size/2+1:seg_size,1:seg_size/2],2)=2*pi*real(ifft(fyxz));  % Partial cumulant density, PBMB (8.11)
t([seg_size/2+1:seg_size,1:seg_size/2],3)=rhoyxz;  % rhoyxz, same indexing as q in col 2, JNM (19)

% Estimate variance of cumulant density and rho estimates
var_fac=4*pi*pi/(seg_size*samp_tot);                             % Factor (2pi/S)(2pi/R)
q_var=var_fac*2*sum(fxxz(1:seg_size/2+1).*fyyz(1:seg_size/2+1)); % Based on PBMB (6.10)
rho_var=1/samp_tot;                                              % Based on JIN (2.40)

cl.type=0;                   % Analysis type
cl.seg_size=seg_size;        % T
cl.seg_tot=seg_tot;          % L
cl.seg_tot_var=seg_tot;      % Effective no of segments, L'
cl.samp_tot=seg_tot*seg_size;% R=LT
cl.samp_rate=samp_rate;      % Sampling rate
cl.dt=deltat;                % Delta t
cl.df=deltaf;                % Delta f
cl.f_c95=0.8512*sqrt(1/cl.seg_tot_var);  % 95% Confidence limit for spectral estimates, PBMB (6.2).
                                         % Confidence interval in plot is TWICE this value
cl.ch_c95=1-0.05^(1/(cl.seg_tot_var-1)); % 95% Confidence limit for partial coherence, PBMB (9.2), JNM (38)
cl.q_c95=1.96*sqrt(q_var);               % 95% Confidence limits for partial cumulant, based on PBMB (6.11).
cl.rho_c95=1.96*sqrt(rho_var);           % 95% Confidence limits for rho, JNM (39)
cl.col_R20=10;               % column in f matrix containing zero lag partial coherence component
cl.col_R2p=11;               % column in f matrix containing forward (positive lag) partial coherence component
cl.col_R2n=12;               % column in f matrix containing reverse (negative lag) partial coherence component
cl.col_rho=3;                % column in t matrix containing rho estimate
cl.R2=R2_rho;                % R2 value, R^2_yx|z, JNM (20)
cl.R2_0=R2_rho_0;            % Component of R2 at zero lag, R^2_yx|z;0. Terms in JNM (22)
cl.R2_p=R2_rho_p;            % Component of R2 in forward direction, positive lag component, R^2_yx|z;+
cl.R2_n=R2_rho_n;            % Component of R2 in reverse direction, negative lag component, R^2_yx|z;-
cl.R2_ch=R2_ch;              % R2 estimated from partial coherence between original processes, JNM (5)
cl.R2_chw=R2_chw;            % R2 estimated from partial coherence between MMSE whitened processes, JNM (18)
cl.R2_fprime_0=R2_fprime_0;  % Component of R2 at zero lag estimated from f'
cl.R2_fprime_p=R2_fprime_p;  % Component of R2 at positive lag estimated from f'
cl.R2_fprime_n=R2_fprime_n;  % Component of R2 at negative lag estimated from f'
cl.N2=0;        % N2, No of events analysed from sp1.
cl.N1=0;        % N1, No of events analysed from sp2.
cl.P2=0;        % P2, mean intensity of events from sp1.
cl.P1=0;        % P1, mean intensity of events from sp2.
cl.opt_str='';  % Options string.
cl.what='';     % Field for plot label.

% Display messages
disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),' sec,  Resolution: ',num2str(cl.df),' Hz.']);
disp(['  R2 (0, p, n): ',num2str(cl.R2),' (',num2str([cl.R2_0 cl.R2_p cl.R2_n]),').'])
