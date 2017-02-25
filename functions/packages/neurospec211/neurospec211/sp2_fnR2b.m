function [f,t,cl]=sp2_fnR2b(dx,dy,fxx,fyy,norm_fac,samp_tot)
% function [f,t,cl]=sp2_fnR2b(dx,dy,fxx,fyy,norm_fac,samp_tot)
%
% Function to implement R2 and directionality measures for bivariate data
% see NOTE below about calling this routine.
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
%
% Inputs
%      dx,dy      2-D matrices of DFT of segments
%      fxx,fyy    Corresponding spectral estimates - used to determine MMSE whitening filters
%      norm_fac   Normalisation factor in spectral estimates, applied as multiplier
%      samp_tot   No of samples in data set.
%
% Outputs
%       f, t, cl  frequency domain matrix, time domain vector, cl structure
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%     Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Halliday DM (2015) Nonparametric directionality measures for time series and point process data,
%     Journal of Integrative Neuroscience, 14(2), In Press. DOI: 10.1142/S0219635215300127
%
% References 1 and 3 referred to in comments as PBMB and JIN, respectively
%
% function [f,t,cl]=sp2_fnR2b(dx,dy,fxx,fyy,norm_fac,samp_tot)
%
% NOTE: This routine is not intended to support analysis of raw data.
% It is intended as a support routine for the 2 channel spectral analysis functions:
% sp2a2_R2_mt.m and sp2_fn2_R2a.m Refer to these functions for further details.

% dx, dy have T rows and L columns
[seg_size,seg_tot]=size(dx);

%Pre-whitening stage - generate MMSE filters, weight for zero frequency set to zero
wx=[0; 1./sqrt(fxx(2:seg_size))];  % JIN (2.23)
wy=[0; 1./sqrt(fyy(2:seg_size))];  % JIN (2.24)

% Apply MMSE pre-whitening filters
dxw=dx;
dyw=dy;
for ind=1:seg_tot
  dxw(:,ind)=dxw(:,ind).*wx;  % JIN (2.25)
  dyw(:,ind)=dyw(:,ind).*wy;  % JIN (2.26)
end

% Re-calculate spectra using pre-whitened processes
fxxw=norm_fac*sum(abs(dxw.*dxw),2);   % Spectra of MMSE whitened process x. This is 1 at all freqs
fyyw=norm_fac*sum(abs(dyw.*dyw),2);   % Spectra of MMSE whitened process y. This is 1 at all freqs
fyxw=norm_fac*sum(dyw.*conj(dxw),2);  % Complex cross spectrum between whitened processes
chyxw=abs(fyxw.*fyxw);                % Coherence from fyxw only, JIN (2.27)
rhoyx=real(ifft(fyxw));               % rhoyx, no factors as 1/T in ifft(), estimate of JIN (2.8)

% Estimate R2 in time domain and separate components:
% zero lag R^2_yx;0, positive lag R^2_yx;+ and negative lag R^2_yx;-
R2_rho=sum(rhoyx.^2);                          % JIN (2.29)
R2_rho_0=rhoyx(1).^2;                          % JIN (2.31) zero lag, bin 1
R2_rho_p=sum(rhoyx(2:seg_size/2).^2);          % JIN (2.32) positive lags, bins 2:(T/2)
R2_rho_n=sum(rhoyx(seg_size/2+1:seg_size).^2); % JIN (2.30) negative lags, bins (T/2)+1:T

% Decompose rho into three components by lag, to estimate f' functions
rhoyx_3=zeros(length(rhoyx),3);
rhoyx_3(1,1)=rhoyx(1);                                         % zero lag
rhoyx_3(2:seg_size/2,2)=rhoyx(2:seg_size/2);                   % positive lags
rhoyx_3(seg_size/2+1:seg_size,3)=rhoyx(seg_size/2+1:seg_size); % negative lags

% Switch back to frequency domain, these are f' estimates
f_prime=fft(rhoyx_3);      % JIN (2.34, 2.35, 2.33):  f'yx;0, f'yx;+, f'yx;-
f_prime2=abs(f_prime).^2;  % These are |f'|^2 values: |f'yx;0|^2, |f'yx;+|^2, |f'yx;-|^2

% Scale to get relative contributions from 0,+,- components.
f_prime2_sum=sum(f_prime2,2);
R2_weight=f_prime2;
R2_weight(:,1)=R2_weight(:,1)./f_prime2_sum;  % Weighting factor in JIN (2.19), range [0 1]
R2_weight(:,2)=R2_weight(:,2)./f_prime2_sum;  % Weighting factor in JIN (2.20), range [0 1]
R2_weight(:,3)=R2_weight(:,3)./f_prime2_sum;  % Weighting factor in JIN (2.18), range [0 1]

%-----------------------------------------------------------------------
% Output matrix f. This is reduced version containing only directionality measures
f=zeros(seg_size/2,7);
f_index=(2:seg_size/2+1)';    % Indexing for output, DC not output.
f(:,1)=fxxw(f_index);         % Column 1 - MMSE whitened spectra process x, JIN (2.4)     
f(:,2)=fyyw(f_index);         % Column 2 - MMSE whitened spectra process y, JIN (2.4)     
f(:,3)=chyxw(f_index);        % Column 3 - Coherence between whitened processes, JIN (2.5)/(2.27)
f(:,4)=angle(fyxw(f_index));  % Column 4 - Phase between whitened processes               
f(:,5)=R2_weight(f_index,1).*chyxw(f_index);   % Column 5 - Zero lag coherence component, JIN (2.19)
f(:,6)=R2_weight(f_index,2).*chyxw(f_index);   % Column 6 - Forward  coherence component, JIN (2.20)
f(:,7)=R2_weight(f_index,3).*chyxw(f_index);   % Column 7 - Reverse  coherence component, JIN (2.18)

% Construct time domain matrix t. Contains only rhoyx estimate, lag  range (-T/2)*dt to (T/2-1)*dt
t([seg_size/2+1:seg_size,1:seg_size/2],1)=rhoyx;  % rhoyx, JIN (2.8)

% Variance of rho
rho_var=1/samp_tot;  % JIN (2.40)

% R2 metrics and confidence limits in structure cl
cl.rho_c95=1.96*sqrt(rho_var); % 95% Confidence limits for rho, JIN (2.41)
cl.R2=R2_rho;                  % R2 value, R^2_yx, JIN (2.29)
cl.R2_0=R2_rho_0;              % Component of R2 at zero lag, R^2_yx;0                                 
cl.R2_p=R2_rho_p;              % Component of R2 in forward direction, positive lag component, R^2_yx;+
cl.R2_n=R2_rho_n;              % Component of R2 in reverse direction, negative lag component, R^2_yx;-
