function psp_rho1(t,cl,lag_tot,lag_neg,label)
% psp_rho1(t,cl,lag_tot,lag_neg,label)
% function to plot rho function from directional analysis in current figure/subplot window
%
% Copyright (C) 2008, David M. Halliday.
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
%     Contact:  contact@neurospec.org
%  
% t,cl     Output from spectral analysis routine.
% lag_tot  Total lag range for time domain including -ve lags (ms).
% lag_neg  Negative lag range for time domain (ms).
% label    Optional title instead of cl.what.
%
%  psp_rho1(t,cl,lag_tot,lag_neg,label)

if (nargin<3)
  error('Not enough input arguments')
end  

if (nargin<4)
  lag_neg=lag_tot/2;
end  
bin_start=cl.seg_size/2+1-round(lag_neg/cl.dt);
bin_stop=bin_start+round(lag_tot/cl.dt)-1;
%Check lag range
[x,y]=size(t);
if (bin_start<1 | bin_start>=bin_stop | bin_stop>x)
  error('Error in requested lag range.');
end

q_zl=[-lag_neg,0;lag_tot-lag_neg,0];
q_ucl=[-lag_neg,cl.rho_c95;lag_tot-lag_neg,cl.rho_c95];
q_lcl=[-lag_neg,-cl.rho_c95;lag_tot-lag_neg,-cl.rho_c95];
plot(t(bin_start:bin_stop,1),t(bin_start:bin_stop,cl.col_rho),'k-',q_zl(:,1),q_zl(:,2),'k--',q_ucl(:,1),q_ucl(:,2),'k-',q_lcl(:,1),q_lcl(:,2),'k-');
axis([-lag_neg,lag_tot-lag_neg,-Inf,Inf]);
xlabel('Lag (ms)')
if (nargin>4)
  title(label);
else
 title(['rho: ',cl.what]);
end  