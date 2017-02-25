% script bar analysis 160105

source_dir = uigetdir('D:\# Projects (Noam)\# SLITE\# DATA\');
D = dir(fullfile(source_dir, '*.txt'));
cd(source_dir);

stimFlag = 0;
cell = {};