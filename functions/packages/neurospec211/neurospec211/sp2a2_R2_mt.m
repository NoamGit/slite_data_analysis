function [f,t,cl,sc] = sp2a2_R2_mt(dat1,dat2,samp_rate,seg_pwr,opt_str)
% [f,t,cl,sc] = sp2a2_R2_mt(dat1,dat2,samp_rate,seg_pwr,opt_str)
%
% Two channel spectral analysis, includes R2 directionality analysis
%  with options including multitaper analysis.
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
%    Contact:  contact@neurospec.org
%
% Input arguments
%  dat1       Channel 1  (input) time series vector.
%  dat2       Channel 2 (output) time series vector.
%  samp_rate  Sampling rate (samples/sec).
%  seg_pwr    Segment length - specified as power of 2.
%  opt_str    Options string.
%
% Output arguments
%  f    column matrix with frequency domain parameter estimates.
%  t    column matrix with      time domain parameter estimates.
%  cl   single structure with scalar values related to analysis.
%  sc   optional matrix for spectral coefficients.
%
% Input Options
%  n  Normalise to unit variance within each segment.
%  r  Rectification - this options requires an argument, valid arguments are:[0, 1, 2].
%                      r0:  Rectify Input  channel  (dat1)
%                      r1:  Rectify Output channel  (dat2)
%                      r2:  Rectify Both   channels (dat1 & dat2)
%  t  linear De-trend - this options requires an argument, valid arguments are:[0, 1, 2].
%                      t0:  De-trend Input  channel  (dat1)
%                      t1:  De-trend Output channel  (dat2)
%                      t2:  De-trend Both   channels (dat1 & dat2)
%  w  Use data window  - this options requires an argument, valid arguments are:[1, 2, 3, 4].
%                      w1: Split cosine taper applied to each segment,  5% tapered at each end.
%                      w2: Split cosine taper applied to each segment, 10% tapered at each end.
%                      w3: Split cosine taper applied to each segment, 25% tapered at each end.
%                      w4: Full  cosine taper applied to each segment, 50% tapered at each end.
%  M  Use multi-taper analysis - requires an argument argument [1, 1.5, 2, 2.5, ..., 5.5, 6].
%                      Option specified as 'Mnw', where nw is the duration x half-bandwidth product, NW.
%                      Number of tapers: K = round(2*nw)-1.
%                      This option is NOT comptabile with 'w' option.
%
% Output parameters
%  f column 1       frequency in Hz.
%  f column 2       Log input/dat1  spectrum.
%  f column 3       Log output/dat2 spectrum.
%  f column 4       Coherence.
%  f column 5       Phase.
%  f column 6 - 12  Frequency domain directionality
%
%  t column 1       Lag in ms.
%  t column 2       Cumulant density estimate
%  t column 3       Time domain directionality, rho estimate
%
%  cl.rho_c95       95% confidence limit for rho estimate
%  cl.R2            R2 value, R^2_yx
%  cl.R2_0          Component of R2 at zero lag, R^2_yx;0
%  cl.R2_p          Component of R2 in forward direction, positive lag component, R^2_yx;+
%  cl.R2_n          Component of R2 in reverse direction, negative lag component, R^2_yx;-
%  cl.col_R20       column in f matrix containing zero lag coherence component
%  cl.col_R2p       column in f matrix containing forward (positive lag) coherence component
%  cl.col_R2n       column in f matrix containing reverse (negative lag) coherence component
%  cl.col_rho       column in t matrix containing rho estimate
%  cl.R2_ch         R2 estimated from coherence between original processes
%  cl.type          Analysis type 0, for compatability with NeuroSpec2.0
%  cl.seg_size      Segment length, T
%  cl.seg_tot       Number of segments, L
%  cl.seg_tot_var   Effective no of segments, used to calculate confidence limits
%  cl.samp_tot      Number of samples analysed, R=LT
%  cl.samp_rate     Sampling rate of data (samples/sec)
%  cl.dt            Time domain bin width (ms)
%  cl.df            Frequency domain bin width (Hz)
%  cl.f_c95         95% confidence limit for Log spectral estimates
%  cl.ch_c95        95% confidence limit for coherence estimate
%  cl.q_c95         95% confidence limit for cumulant density estimate
%  cl.MT            Flag indicating multi-taper (MT) estiamte - 0: No, 1: Yes
%  cl.MT_NW         NW for MT estimates
%  cl.MT_tap        No of tapers in MT estiamte
%  cl.MT_eig        Eigenvalues for individual eigen spectra in MT estimate
%  cl.N1, cl.N2     Spike counts for point process data
%  cl.P1, cl.P2     Spike rates for point process data
%  cl.opt_str       Copy of options string
%  cl.what          Text label, used in plotting routines
%
%  sc column 1      f11.
%  sc column 2      f22.
%  sc column 3      f21 (Complex).
%
% References:
% 1. Halliday D.M., Rosenberg J.R., Amjad A.M., Breeze P., Conway B.A. & Farmer S.F.
%    Progress in Biophysics and molecular Biology, 64, 237-278, 1995.
% 2. Bloomfield, P. Fourier Analysis of Time Series: An Introduction.
%    2nd edition. Wiley, New York, 2000.
% 3. Halliday DM, Journal of Integrative Neuroscience, 14(2), In Press. DOI: 10.1142/S0219635215300127
%
% [f,t,cl,sc] =  sp2a2_R2_mt(dat1,dat2,samp_rate,seg_pwr,opt_str)

% Check numbers of arguments
if (nargin<4)
  error(' Not enough input arguments');
end
if (nargout<3)
  error(' Not enough output arguments');
end

% Check for single column data
[nrow,ncol]=size(dat1);
if (ncol~=1)
  error(' Input NOT single column: dat1')
end
[nrow,ncol]=size(dat2);
if (ncol~=1)
  error(' Input NOT single column: dat2')
end

pts_tot=length(dat1);           % Determine size of data vector.
if (length(dat2)~=pts_tot)      % Check that input vectors are equal length. 
  error (' Unequal length data arrays');
end

if (max(size(samp_rate)) ~= 1)
  error(' Type 0 - Non scalar value for: samp_rate');
end
if (max(size(seg_pwr)) ~= 1)
  error(' Type 0 - Non scalar value for: seg_pwr');
end
seg_size=2^seg_pwr;                    % DFT segment length (T).
seg_tot=fix(pts_tot/seg_size);         % Number of complete segments (L).
samp_tot=seg_tot*seg_size;             % Number of samples to analyse: R=LT.
seg_samp_start=(1:seg_size:samp_tot)'; % Start times for each segment.
seg_samp_no(1:seg_tot,1)=seg_size;     % Fixed number of data points per segment, T.

% Create data matrices, T rows, L columns.
rd1=zeros(seg_size,seg_tot);
rd2=zeros(seg_size,seg_tot);

%Maximum number of tapers & seg_pwr limits for Multi-taper analysis
MT_tap_max=11;
MT_seg_pwr_min=6;
MT_seg_pwr_max=14;

% Set option defauts
flags.sp_type=0;  % Type 0 analysis
flags.wind=0;     % Set defaults - options off.
flags.MT=0;
norm_chan=0;
trend_chan_1=0; 
trend_chan_2=0;
rect_chan_1=0;
rect_chan_2=0;
if (nargin<5)
  opt_str='';
end
options=deblank(opt_str);
while (any(options))              % Parse individual options from string.
  [opt,options]=strtok(options);
  optarg=opt(2:length(opt));      % Determine option argument.
  switch (opt(1))
    case 'n'             % Normalisation to unit variance within each segment.
      norm_chan=1;
    case 'r'             % Rectification option.        
      i=str2num(optarg);
    if (i<0 | i>2)
      error(['error in option argument -- r',optarg]);
    end  
    if (i~=1)
      rect_chan_1=1;     % Rectify ch 1.
    end  
    if (i>=1)
      rect_chan_2=1;     % Rectify ch 2.
    end  
    case 't'             % Linear de-trend option.
      i=str2num(optarg);
    if (i<0 | i>2)
      error(['error in option argument -- t',optarg]);
    end  
    if (i~=1)
      trend_chan_1=1;    % De-trend ch 1.
    end  
    if (i>=1)
      trend_chan_2=1;    % De-trend ch 2.
    end  
    case 'w'             % Data window (taper) option
      i=str2num(optarg);
    if (i<1 | i>4)
      error(['error in option argument -- w',optarg]);
    end
    flags.wind=i;
    case 'M'             % Multi-taper option
      if (seg_pwr<MT_seg_pwr_min | seg_pwr>MT_seg_pwr_max)
        error(['seg_pwr outside allowed range for Multi-taper otpion:',num2str([MT_seg_pwr_min MT_seg_pwr_max])]);
      end
      MT_flag=1;
      NW=str2num(optarg);
      i=round(2*NW)-1; % K ~ 2NW - 1
      if (i<1 | i>MT_tap_max)
        error(['error in option argument -- M',optarg]);
      end
      if MT_flag
        flags.MT=1;
        flags.seg_pwr=seg_pwr;
        flags.NW=NW;
        flags.MT_tap=i;
      end 
    otherwise
      error (['Illegal option -- ',opt]);  % Illegal option.
  end
end

% Check MT option compatability
if (flags.MT & flags.wind)
  disp('Warning - Data window option not compatible with MT analysis. Option disabled')
  flags.wind=0;
end

if (trend_chan_1 | trend_chan_2)       % Index for fitting data with polynomial.
  trend_x=(1:seg_size)';
end

for ind=1:seg_tot                     % Loop across columns/segments.
  seg_pts=seg_samp_no(ind);           % No of data points in segment.
  seg_start=seg_samp_start(ind);      % Start sample in data vector.
  seg_stop=seg_start+seg_pts-1;       % Stop  sample in data vector.
  dat_seg1=dat1(seg_start:seg_stop);  % Extract segment from dat1.
  dat_seg2=dat2(seg_start:seg_stop);  % Extract segment from dat2.
  md1=mean(dat_seg1);                 % Mean of segment from dat1.
  md2=mean(dat_seg2);                 % Mean of segment from dat2.

  rd1(1:seg_pts,ind)=dat_seg1-md1;    % Subtract mean from ch 1.
  rd2(1:seg_pts,ind)=dat_seg2-md2;    % Subtract mean from ch 2.
  if rect_chan_1
    rd1(1:seg_pts,ind)=abs(rd1(1:seg_pts,ind));            % Rectification of ch 1 (Full wave).
  end
  if rect_chan_2
    rd2(1:seg_pts,ind)=abs(rd2(1:seg_pts,ind));            % Rectification of ch 2 (Full wave).  
  end
  if trend_chan_1                               % Linear trend removal.
    p=polyfit(trend_x(1:seg_pts,1),rd1(1:seg_pts,ind),1);                 % Fit 1st order polynomial.
    rd1(1:seg_pts,ind)=rd1(1:seg_pts,ind)-p(1)*trend_x(1:seg_pts,1)-p(2); % Subtract from ch 1.
  end  
  if trend_chan_2                               % Linear trend removal.
    p=polyfit(trend_x(1:seg_pts,1),rd2(1:seg_pts,ind),1);                 % Fit 1st order polynomial.
    rd2(1:seg_pts,ind)=rd2(1:seg_pts,ind)-p(1)*trend_x(1:seg_pts,1)-p(2); % Subtract from ch 2.
  end
  if norm_chan
    rd1(1:seg_pts,ind)=rd1(1:seg_pts,ind)/std(rd1(1:seg_pts,ind));
    rd2(1:seg_pts,ind)=rd2(1:seg_pts,ind)/std(rd2(1:seg_pts,ind));
  end
end

% Call sp2_fn2_R2a MT based spectral estimation routine, with R2 directionality analysis
if (nargout>3)
  [f,t,cl,sc]=sp2_fn2_R2a(rd1,rd2,seg_samp_no,samp_rate,flags);
else  
  [f,t,cl]=sp2_fn2_R2a(rd1,rd2,seg_samp_no,samp_rate,flags);
end 

% Set additional data and type dependent elements in cl structure.
cl.N1=0;            % N1, No of events in ch 1. (zero for TS data)
cl.N2=0;            % N2, No of events in ch 2.          "
cl.P1=0;            % P1, mean intensity ch1.            "
cl.P2=0;            % P2, mean intensity ch2.            "
cl.opt_str=opt_str; % Copy of options string.
cl.what='';         % Field for plot label.

% Display No of segments & resolution
disp(['Segments: ',num2str(seg_tot),', Segment length: ',num2str(seg_size/samp_rate),' sec,  Resolution: ',num2str(cl.df),' Hz.']);
disp(['  R2 (0, p, n): ',num2str(cl.R2),' (',num2str([cl.R2_0 cl.R2_p cl.R2_n]),').'])
