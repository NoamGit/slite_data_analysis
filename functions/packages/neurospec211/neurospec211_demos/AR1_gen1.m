function [dat]=AR1_gen1(pts_tot)
% function [dat]=AR1_gen1(pts_tot)
% 
% Function to generate data for testing R2 directionality analysis, 
% generates a single AR(1) process, with pts_tot samples.
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

% Coefficent in AR(1) model.
% Specified as no of steps per time constant.
% Sensible range would be 5-25 steps.
tau_steps=5;
% AR(1) coefficient a.
a=exp(-1/tau_steps);

% Ouput process has desired variance of 1
% Calculate variance of input noise sequence, w
var_w=1-a^2;
std_w=sqrt(var_w);

% Generate signal noise sequence, w(k)
w=std_w*randn(pts_tot+1,1);

% Vector for signal, x
x=zeros(pts_tot+1,1);

% Generate x recursively as AR(1) process
for k=2:pts_tot+1
  x(k)=a*x(k-1)+w(k-1);  % AR(1) equation
end

% Output from k=2 onwards
dat=x(2:pts_tot+1);
