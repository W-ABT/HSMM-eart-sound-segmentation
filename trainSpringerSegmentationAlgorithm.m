% function [logistic_regression_B_matrix, pi_vector, total_obs_distribution] = trainSpringerSegmentationAlgorithm(PCGCellArray, annotationsArray, Fs, figures)
%
% Training the Springer HMM segmentation algorithm. Developed for use in
% the paper:
%训练Springer-HMM分割算法。为在论文中使用而开发：
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% Inputs:
% PCGCellArray: A 1XN cell array of the N audio signals. For evaluation
% purposes, these signals should be from a distinct training set of
% recordings, while the algorithm should be evaluated on a separate test
% set of recordings, which are recorded from a completely different set of
% patients (for example, if there are numerous recordings from each
% patient).
%PCGCellArray:由n个音频信号组成的1xn单元阵列。出于评估目的，这些信号应来自不同的记录训练集，
%而算法应在记录的单独测试集上进行评估，记录来自完全不同的患者集（例如，如果每个患者有许多记录）。
% annotationsArray: a Nx2 cell array: position (n,1) = the positions of the
% R-peaks and postion (n,2) = the positions of the end-T-waves
% (both in SAMPLES)
%annotationsArray:Nx2单元阵列：位置（n，1）=R峰的位置，位置（n，2）=末端T波的位置（均在样本中）
% Fs: The sampling frequency of the PCG signals
%FS:PCG信号的采样频率
% figures (optional): boolean variable dictating the disaplay of figures.
%图形（可选）：指示图形不显示的布尔变量。
%% Outputs:
% logistic_regression_B_matrix:逻辑回归矩阵
% pi_vector:π矢量
% total_obs_distribution:总体OBS分布
% As Springer et al's algorithm is a duration dependant HMM, there is no
% need to calculate the A_matrix, as the transition between states is only
% dependant on the state durations.
%由于springer等人的算法是一个依赖于持续时间的hmm，因此不需要计算a_矩阵，
%因为状态之间的转换仅依赖于状态持续时间。
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


function [logistic_regression_B_matrix, pi_vector, total_obs_distribution] = trainSpringerSegmentationAlgorithm(PCGCellArray, annotationsArray, Fs, figures)

%% Options

if(nargin < 4)
    figures = false;
end


numberOfStates = 4;%状态数量有4个，设隐藏态定义为ζ=[ζ1，ζ2，…，ζN]，
%其中N是状态总数。其中N=4，ζ1为s1，ζ2为收缩期，ζ3为s2，ζ4为舒张期。
numPCGs = length(PCGCellArray);%

% A matrix of the values from each state in each of the PCG recordings:
%每个PCG记录中每个状态的值矩阵：
state_observation_values = cell(numPCGs,numberOfStates);
%cell：返回一个个a*b的数组，状态观测值

for PCGi = 1:length(PCGCellArray)
    PCG_audio = PCGCellArray{PCGi};
    
    S1_locations = annotationsArray{PCGi,1};%S1定位
    S2_locations = annotationsArray{PCGi,2};%S2定位
    
    [PCG_Features, featuresFs] = getSpringerPCGFeatures(PCG_audio, Fs);
    
    PCG_states = labelPCGStates(PCG_Features(:,1),S1_locations, S2_locations, featuresFs);
    
    
    %% Plotting assigned states:
    %绘制指定状态：
    if(figures)
        figure('Name','Assigned states to PCG');
        
        t1 = (1:length(PCG_audio))./Fs;
        t2 = (1:length(PCG_Features))./featuresFs;
        
        plot(t1, PCG_audio, 'k-');
        hold on;
        plot(t2, PCG_Features, 'b-');
        plot(t2, PCG_states, 'r-');
        
        legend('Audio','Features','States');
        pause();
    end
    
    
    
    %% Group together all observations from the same state in the PCG recordings:
    %将PCG记录中来自同一状态的所有观察结果组合在一起：
    for state_i = 1:numberOfStates
        state_observation_values{PCGi,state_i} = PCG_Features(PCG_states == state_i,:);
    end
end

% Save the state observation values to the main workspace of Matlab for
% later investigation if needed:如果需要，将状态观测值保存到Matlab的主工作空间中，以便以后进行研究：
assignin('base', 'state_observation_values', state_observation_values)

%% Train the B and pi matrices after all the PCG recordings have been labelled:
%在对所有PCG记录进行标记后，训练B和pi矩阵：
[logistic_regression_B_matrix, pi_vector, total_obs_distribution] = trainBandPiMatricesSpringer(state_observation_values);

