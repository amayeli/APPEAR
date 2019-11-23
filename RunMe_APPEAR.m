% Copyright (C) 2019 LaureateInstitute for Brain Research
% For any questtions please contact: amayeli@laureateinstitute.org
%
% Description: 
% Setting Parameters to Run APPEAR
setting_file='SettingsExample.log';
input_dir='Example';
output_dir='Results';
output_file='Out_File';
scan_time=480;
removal_time=0;
TR=2;
slice_per_TR=39;
slice_marker_per_TR=39;
cd Example
EEG_file = dir('*.eeg');
tic;
file = 'AA181-T0-REST-R1-RAW.eeg';
C = strsplit(file,'.');
input_file=C{1};

APPEAR(setting_file,input_dir,input_file,output_dir,output_file,scan_time,removal_time,TR,slice_per_TR,slice_marker_per_TR)
t=toc;