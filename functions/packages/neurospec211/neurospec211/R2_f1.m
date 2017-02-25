function [cl_out,R2_out]=R2_f1(f,cl_in,max_freq)
% function [cl_out,R2_out]=R2_f1(f,cl_in,max_freq)
%
%  Calculation of R2 directionality metrics over restricted frequency range
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
% Input arguments
%  f         Frequency matrix from sp2a2_R2 or sp2a2_R2_mt analysis
%  cl_in     cl structure     from sp2a2_R2 or sp2a2_R2_mt analysis
%  max_freq  Frequency to use in calculating R2 directionality metrics
%
%
% Output arguments (Optional)
%  cl_out   Structure with R2 directionality metrics over reduced frequency range
%  R2_out   Vector of R2 directionality metrics over reduced frequency range
%
%  cl_out.R2_ch       R2 value, R^2_yx, over full frequency range
%  cl_out.R2_ch_f     R2 value, R^2_yx, over reduced frequency range
%  cl_out.R2_0_f      Component of R2 at zero lag over reduced frequency range, R^2_yx;0,fa
%  cl_out.R2_p_f      Component of R2 in forward direction over reduced frequency range, R^2_yx;+,fa
%  cl_out.R2_n_f      Component of R2 in reverse direction over reduced frequency range, R^2_yx;-,fa
%  cl_out.R2_per_0    Percentage of R2 due to zero lag component, over reduced frequency range
%  cl_out.R2_per_p    Percentage of R2 in forward direction, over reduced frequency range
%  cl_out.R2_per_n    Percentage of R2 in reverse direction, over reduced frequency range
%  cl_out.R2_freq     Reduced frequency range used to calculate R2 metrics, fa
%  cl_out.R2_freq_pts Number of frequency bins used in calculation of R2 metrics
%  cl_out.seg_size    Segment length, T
%  cl_out.samp_rate   Sampling rate, samples/sec
%  cl_out.df          Frequency resolution, delta f
% 
% Reference:
% 1. Halliday DM (2015) Nonparametric directionality measures for time series and point process data,
%     Journal of Integrative Neuroscience, 14(2), In Press.   DOI: 10.1142/S0219635215300127
%
% function [cl_out,R2_out]=R2_f1(f,cl_in,max_freq)


% Restricted frequency range decomposition using parameter freq_range
freq_pts=round(max_freq/cl_in.df);

% R2 decomposition over restricted frequency range - includes contribution of zero from zero frequency
R2_ch_f=(2/(2*freq_pts+1))*sum(f(1:freq_pts,4));  % JIN (2.36)
R2_0_f=(2/(2*freq_pts+1))*sum(f(1:freq_pts,10));  % JIN (2.37)
R2_p_f=(2/(2*freq_pts+1))*sum(f(1:freq_pts,11));  % JIN (2.37)
R2_n_f=(2/(2*freq_pts+1))*sum(f(1:freq_pts,12));  % JIN (2.37)

% Percentage R2 in zero lag, forward, reverse direction
R2_per_0=100*R2_0_f/R2_ch_f;  % Percentage at zero lag
R2_per_p=100*R2_p_f/R2_ch_f;  % Percentage in forward direction
R2_per_n=100*R2_n_f/R2_ch_f;  % Percentage in reverse direction

% Display messages
disp(['R2_f1 - Freq: ',num2str(max_freq),',  Pts: ',num2str(freq_pts),'.  R2: ',num2str(R2_ch_f),', (0, p, n): ',num2str([R2_0_f, R2_p_f, R2_n_f]),'.  Sum: ',num2str([R2_0_f+R2_p_f+R2_n_f])])

if (nargout>0)
  cl_out.R2_ch=cl_in.R2_ch;
  cl_out.R2_ch_f=R2_ch_f;
  cl_out.R2_0_f=R2_0_f;
  cl_out.R2_p_f=R2_p_f;
  cl_out.R2_n_f=R2_n_f;
  cl_out.R2_per_0=R2_per_0;
  cl_out.R2_per_p=R2_per_p;
  cl_out.R2_per_n=R2_per_n;

  cl_out.R2_freq=max_freq;
  cl_out.R2_freq_pts=freq_pts;
  cl_out.seg_size=cl_in.seg_size;
  cl_out.samp_rate=cl_in.samp_rate;
  cl_out.df=cl_in.df;
end

if (nargout>1)
  R2_out=[R2_ch_f; R2_0_f; R2_p_f; R2_n_f];
end 
