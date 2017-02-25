function [dat]=AR4_gen1(pts_tot)
% function [dat]=AR4_gen1(pts_tot)
% 
% Function to generate data for testing R2 directionality analysis
% Based on univariate AR(4) example from Percival & Walden, 1993, Equation 46a.
% Two univariate AR(4) processes generated with connection from 1->2 only.
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
% Reference for AR(4) equations:
%  Percival DB & Walden AT (1993) Spectral Analysis for Physical Applications,
%  Cambridge University Press, Page 46.

AR_order=4;  % Order of AR processes
% Coefficients from equation 46a in Percival & Walden, 1993
a4=[2.7607 -3.8106 2.6535 -0.9238];

% Extend length to include AR order
pts_tot=pts_tot+AR_order;
x1=zeros(pts_tot,1);
x2=zeros(pts_tot,1);

% Separate N(0,1) white noise processes
w1=randn(pts_tot,1);
w2=randn(pts_tot,1);

% Weight to use for coupling AR processes
AR_w=0.01;

% Generate two AR(4) signals with coupling only from 1->2
for n=AR_order+1:pts_tot
  x1(n)=a4*x1(n-1:-1:n-4)+w1(n);              % Process 1
  x2(n)=a4*x2(n-1:-1:n-4)+AR_w*x1(n-1)+w2(n); % Process 2, with input from Process 1
end

% Standardise mean/variance
x1=(x1-mean(x1))/std(x1);
x2=(x2-mean(x2))/std(x2);

% Output samples after initial AR_order values
dat=[x1(AR_order+1:pts_tot) x2(AR_order+1:pts_tot)];
