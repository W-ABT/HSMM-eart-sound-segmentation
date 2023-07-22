%% Example Springer script
% A script to demonstrate the use of the Springer segmentation algorithm

%% Copyright (C) 2016  David Springer
% dave.springer@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%
clear;
clc;

%% Load the default options:
% These options control options such as the original sampling frequency of
% the data, the sampling frequency for the derived features and whether the
% mex code should be used for the Viterbi decoding:
%����Ĭ��ѡ�
%��Щѡ������������ݵ�ԭʼ����Ƶ�ʡ����������Ĳ���Ƶ���Լ��Ƿ�Ӧ��mex������viterbi�����ѡ�
springer_options = default_Springer_HSMM_options;

%% Load the audio data and the annotations:
% These are 6 example PCG recordings, downsampled to 1000 Hz, with
% annotations of the R-peak and end-T-wave positions.
%������Ƶ���ݺ���ע��
%����6��pcg��¼���²�����1000hz������r���t��ĩ��λ�õ�ע�͡�
load('example_data.mat');
%% Split the data into train and test sets:
% Select the 5 recordings for training and a sixth for testing:
%�����ݷ�Ϊѵ�����Ͳ��Լ���
%ѡ��5����¼������ѵ����6��¼�����в��ԣ�
training_indices = [1, 47, 361, 402, 572];
train_recordings = example_data.example_audio_data(training_indices);
train_annotations = example_data.example_annotations(training_indices,:);

% prompt = 'x\n';
% x = input(prompt);
test_index = [1];

test_recordings = example_data.example_audio_data(test_index);
test_annotations = example_data.example_annotations(test_index,:);

%% Train the HMM:ѵ��HMM
[B_matrix, pi_vector, total_obs_distribution] = trainSpringerSegmentationAlgorithm(train_recordings,train_annotations,springer_options.audio_Fs, false);

%% Run the HMM on an unseen test recording:
% And display the resulting segmentation
%�ڿ������Ĳ��Լ�¼������hmm��
%����ʾ����ֶ�
numPCGs = length(test_recordings);

for PCGi = 1:numPCGs
    [assigned_states,PCG_states,PCG_Features] = runSpringerSegmentationAlgorithm(test_recordings{PCGi}, springer_options.audio_Fs, B_matrix, pi_vector, total_obs_distribution, test_annotations,true);
end

accuracy = length(find(assigned_states == PCG_states))/length(PCG_states)*100;