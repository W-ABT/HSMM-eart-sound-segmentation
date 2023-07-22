% function assigned_states = runSpringerSegmentationAlgorithm(audio_data, Fs, B_matrix, pi_vector, total_observation_distribution, figures)
%
% A function to assign states to a PCG recording using a duration dependant
% logisitic regression-based HMM, using the trained B_matrix and pi_vector
% trained in "trainSpringerSegmentationAlgorithm.m". Developed for use in
% the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%一种函数，使用基于持续时间的逻辑回归HMM，
%使用“trainSpringerSegmentationAlgorithm.m”中训练的B_矩阵和pi_向量，将状态分配给PCG记录。
%% INPUTS:
% audio_data: The audio data from the PCG recording
% Fs: the sampling frequency of the audio recording
% B_matrix: the observation matrix for the HMM, trained in the 
% "trainSpringerSegmentationAlgorithm.m" function
% pi_vector: the initial state distribution, also trained in the 
% "trainSpringerSegmentationAlgorithm.m" function
% total_observation_distribution, the observation probabilities of all the
% data, again, trained in trainSpringerSegmentationAlgorithm.
% figures: (optional) boolean variable for displaying figures
%
%% OUTPUTS:
% assigned_states: the array of state values assigned to the original
% audio_data (in the original sampling frequency).
%分配状态：分配给原始音频数据的状态值数组（以原始采样频率）。
%
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

function [assigned_states,PCG_states,PCG_Features] = runSpringerSegmentationAlgorithm(audio_data, Fs, B_matrix, pi_vector, total_observation_distribution, annotationsArray,figures)


if(nargin < 6)
    figures = false;
end

%% Get PCG Features:

[PCG_Features, featuresFs] = getSpringerPCGFeatures(audio_data, Fs);

%% Get PCG heart rate

[heartRate, systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs);

[~, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, pi_vector, B_matrix, total_observation_distribution, heartRate, systolicTimeInterval, featuresFs);

assigned_states = expand_qt(qt, featuresFs, Fs, length(audio_data));

%% 真实状态
[PCG1_Features, features1Fs] = Copy_of_getSpringerPCGFeatures(audio_data, Fs);
    S1_locations = annotationsArray{1,1};%S1定位
    S2_locations = annotationsArray{1,2};%S2定位

    states = labelPCGStates(PCG1_Features(:,1),S1_locations, S2_locations, features1Fs);
    PCG_states = expand_qt(states, features1Fs, Fs, length(audio_data));
% 出图，呈现出切分后图像与原波形图像经行对比

 %if(figures)
 %   figure('Name','Derived state sequence');
 %   t1 = (1:length(audio_data))./Fs;
 %   plot(t1,normalise_signal(audio_data),'k');
   % if(a)
   % hold on;
   % else
   %    hold off;
   % end
 %  hold on;
 %   plot(t1,assigned_states,'r--');
 %   xlabel('Time (s)');
 %   legend('Audio data', 'Derived states');
 %   grid minor;
end




