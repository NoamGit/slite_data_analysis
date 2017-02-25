function psp_ch1(f,cl,freq,ch_max,label)
% psp_ch1(f,cl,freq,ch_max,label)
% function to plot coherence in current figure/subplot window
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
% f,cl     Output from spectral analysis routine.
% freq     Frequency limit for plotting (Hz).
% ch_max   Maximum value of y axis (Optional).
% label    Optional title instead of cl.what.
%
% psp_ch1(f,cl,freq,ch_max,label)

if (nargin<3)
  error('Not enough input arguments')
end  

%freq_pts=round(freq/cl.df);
%Updated to support type 2 analysis
f_start=f(1,1,1);
freq_pts=round((freq-f_start)/cl.df)+1;

%Check freq range
[x,y]=size(f);
if (freq_pts>x)
  error('Requested frequency range too large.');
end
if (freq_pts<1)
  error('Requested frequency range too narrow.');
end

f_max=f(freq_pts,1);
ch_cl=[f_start,cl.ch_c95;f_max,cl.ch_c95];
plot(f(1:freq_pts,1),f(1:freq_pts,4),'k-',ch_cl(:,1),ch_cl(:,2),'k--');
if (nargin<4)
  axis([0,freq,0,Inf]);
else
  axis([0,freq,0,ch_max]);
end  
xlabel('Frequency (Hz)')
if (nargin>4)
  title(label);
else
 title(['coh: ',cl.what]);
end  