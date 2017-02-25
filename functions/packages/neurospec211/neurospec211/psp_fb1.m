function psp_fb1(f,cl,freq,label)
% psp_fb1(f,cl,freq,label)
% Function to plot log spectra in current figure/subplot window
% Plots spectra in column 3 of f (output signal)
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
% label    Optional title instead of cl.what.
%
% psp_fb1(f,cl,freq,label)

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
if (cl.P2>0)
 % Point process spectrum
  f_cl=[f_start,log10(cl.P2);f_max,log10(cl.P2);
        f_start,log10(cl.P2)+cl.f_c95;f_max,log10(cl.P2)+cl.f_c95;
        f_start,log10(cl.P2)-cl.f_c95;f_max,log10(cl.P2)-cl.f_c95];
  plot(f(1:freq_pts,1),f(1:freq_pts,3),'k-',f_cl(1:2,1),f_cl(1:2,2),'k--',f_cl(3:4,1),f_cl(3:4,2),'k-',f_cl(5:6,1),f_cl(5:6,2),'k-');
  axis([0,freq,-Inf,Inf]);
else
 % Time series spectrum
  f_cl=[f_max*0.95,max(f(1:freq_pts,3));
        f_max*0.95,max(f(1:freq_pts,3))-2*cl.f_c95];
  plot(f(1:freq_pts,1),f(1:freq_pts,3),'k-',f_cl(:,1),f_cl(:,2),'k-');
  axis([0,freq,-Inf,Inf]);
end  

xlabel('Frequency (Hz)')
if (nargin>3)
  title(label);
else
 title(['fb: ',cl.what]);
end  