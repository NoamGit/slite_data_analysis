function psp_ch1_R2(f,cl,freq,ch_max,label)
% psp_ch1_R2(f,cl,freq,ch_max,label)
% function to plot coherence and directional components in current figure/subplot window
%  
% Copyright (C) 2015, David M. Halliday.
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
% psp_ch1_R2(f,cl,freq,ch_max,label)

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

% Ordinary coherence - black line, drawn first
plot(f(1:freq_pts,1),f(1:freq_pts,4),'k-')
hold on
plot(f(:,1),f(:,cl.col_R20),'Color',[0.5 0.5 0.5])  % zero lag components 50% in gray
plot(f(:,1),f(:,cl.col_R2p),'r')                    % forward direction in red
plot(f(:,1),f(:,cl.col_R2n),'b')                    % reverse direction in blue
plot([0 freq],[cl.ch_c95 cl.ch_c95],'k--');         % Confidence limit
hold off
if (nargin<4)
  axis([0,freq,0,Inf]);
else
  axis([0,freq,0,ch_max]);
end  
xlabel('Frequency (Hz)')
if (nargin>4)
  title(label);
else
 title(['coh dir: ',cl.what]);
end
