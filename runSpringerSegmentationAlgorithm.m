% function assigned_states = runSpringerSegmentationAlgorithm(audio_data, Fs, B_matrix, pi_vector, total_observation_distribution, figures)
%
% A function to assign states to a PCG recording using a duration dependant
% logisitic regression-based HMM, using the trained B_matrix and pi_vector
% trained in "trainSpringerSegmentationAlgorithm.m". Developed for use in
% the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%һ�ֺ�����ʹ�û��ڳ���ʱ����߼��ع�HMM��
%ʹ�á�trainSpringerSegmentationAlgorithm.m����ѵ����B_�����pi_��������״̬�����PCG��¼��
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
%����״̬�������ԭʼ��Ƶ���ݵ�״ֵ̬���飨��ԭʼ����Ƶ�ʣ���
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

%% ��ʵ״̬
[PCG1_Features, features1Fs] = Copy_of_getSpringerPCGFeatures(audio_data, Fs);
    S1_locations = annotationsArray{1,1};%S1��λ
    S2_locations = annotationsArray{1,2};%S2��λ

    states = labelPCGStates(PCG1_Features(:,1),S1_locations, S2_locations, features1Fs);
    PCG_states = expand_qt(states, features1Fs, Fs, length(audio_data));
% ��ͼ�����ֳ��зֺ�ͼ����ԭ����ͼ���жԱ�

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




