% NeuroSpec
% Version 2.11  01-Jul-2016
%
%Core routines
%  sp2a2_R2     - Unconditional non-parametric directionality analysis between two signals.
%  sp2a2_R2_mt  - Unconditional non-parametric directionality analysis, with additional options.
%  sp2_fn2_R2a  - Function called by sp2a2_R2_mt, performs spectral analysis.
%  sp2_fnR2b    - Function called by sp2_fn2_R2a, performs directionality analysis.
%  sp2a2_R2_pc1 - Conditional non-parametric directionality analysis using three signals.
%
%R2 metrics over reduced frequency range
%  R2_f1        - Calculation of directionality metrics over reduced frequency range.
%
%Main Plotting Routine
%  psp2_R2      - Plots results from core routines, including directionality in time and frequency domains.
%
%Other Plotting Routines
%  psp_fa1      - Plots input log spectral estimate (from NeuroSpec 2.0).
%  psp_fb1      - Plots output log spectral estimate (from NeuroSpec 2.0).
%  psp_ch1      - Plots coherence estimate (from NeuroSpec 2.0).
%  psp_ph1      - Plots phase estimate (from NeuroSpec 2.0).
%  psp_q1       - Plots cumulant density estimate (from NeuroSpec 2.0).
%  psp_ch1_R2   - Plots coherence estimate with directional components.
%  psp_rho1     - Plots estimate of time domain directionality.
%
%Demonstration scripts
%  R2_cn_demo1       - Demonstrates spike train analysis using simulated cortical neuron data.
%  R2_ts_demo1       - Demonstrates time series analysis on white and non-white data with known directionality.
%  R2_AR4_demo1      - Demonstrates time series analysis on coupled AR(4) signals.
%  R2_pc_table1_demo - Demonstrate conditional analysis with three signals, generates table of estimated values.
%  R2_pc_table2_demo - Demonstrate conditional analysis with four signals, generates table of estimated values.
%  R2_cn_pc_demo1    - Demonstrate conditional analysis using simulated cortical neuron data.
%  R2_ts_pc_demo1    - Demonstrate conditional analysis using mixtures of Gaussian signals with delays.

%Copyright (C) 2016, David M. Halliday.
%This file is part of NeuroSpec.
%
%   NeuroSpec is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   NeuroSpec is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with NeuroSpec; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%   NeuroSpec is available at:  http://www.neurospec.org/
%    Contact:  contact@neurospec.org
