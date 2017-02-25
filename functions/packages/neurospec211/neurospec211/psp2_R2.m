function psp2_R2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)
% psp2_R2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)
% Function to plot output of 2 channel spectral analysis,
%  from type 0, type 1 or single offset type 2.
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
% f,t,cl   Output from spectral analysis routine.
% freq     Frequency limit for plotting (Hz).
% lag_tot  Total lag range for time domain including -ve lags (ms).
% lag_neg  Negative lag range for time domain (ms).    
% ch_max   Upper limit for coherence estimte.
% label    Label for subplots, used in preference to cl.what
%
% This includes R2 measures plotting: R20, R2+, R2-, rho
%
% psp2_R2(f,t,cl,freq,lag_tot,lag_neg,ch_max,label)

% Check numbers of arguments. 
if (nargin<7)
  error('Not enough input arguments');
end

if (nargin<8)
  subplot(4,2,1)
  psp_fa1(f,cl,freq)
  subplot(4,2,2)
  psp_fb1(f,cl,freq)
  subplot(4,2,3)
  psp_ch1(f,cl,freq,ch_max)
  subplot(4,2,4)
  psp_ph1(f,cl,freq)
  subplot(4,2,5)
  psp_ch1_R2(f,cl,freq,ch_max)
	subplot(4,2,6)
  psp_q1(t,cl,lag_tot,lag_neg)
	subplot(4,2,8)
  psp_rho1(t,cl,lag_tot,lag_neg)
else  
  subplot(4,2,1)
  psp_fa1(f,cl,freq,label)
  subplot(4,2,2)
  psp_fb1(f,cl,freq,label)
  subplot(4,2,3)
  psp_ch1(f,cl,freq,ch_max,label)
  subplot(4,2,4)
  psp_ph1(f,cl,freq,label)
  subplot(4,2,5)
  psp_ch1_R2(f,cl,freq,ch_max,label)
	subplot(4,2,6)
  psp_q1(t,cl,lag_tot,lag_neg,label)
	subplot(4,2,8)
  psp_rho1(t,cl,lag_tot,lag_neg,label)
end
