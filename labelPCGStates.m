 function states = labelPCGStates(envelope,s1_positions, s2_positions, samplingFrequency, figures)
%
% This function assigns the state labels to a PCG record. 
% This is based on ECG markers, dervied from the R peak and end-T wave locations.
%
%此函数用于为PCG记录分配状态标签。
%这是基于心电图标记物，从R波峰值和T波末端位置进行分析。
%
%% Inputs:
% envelope: The PCG recording envelope (found in getSchmidtPCGFeatures.m)
% s1_positions: The locations of the R peaks (in samples)
% s2_positions: The locations of the end-T waves (in samples)
% samplingFrequency: The sampling frequency of the PCG recording
% figures (optional): boolean variable dictating the display of figures
% 
% envelope:PCG记录信封（可在getSchmidtPCGFeatures.m中找到）
% s1_位置：R峰的位置（样品中）
% s2_位置：端T波的位置（样本中）
% 采样频率：PCG记录的采样频率
% 图形（可选）：指示图形显示的布尔变量
% 
%% Output:
% states: An array of the state label for each sample in the feature
% vector. The total number of states is 4. Therefore, this is an array of
% values between 1 and 4, such as: [1,1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,1,1,1],
% illustrating the "true" state label for each sample in the features.
% State 1 = S1 sound
% State 2 = systole
% State 3 = S2 sound
% State 4 = diastole
% 
% 状态：特征向量中每个样本的状态标签数组。状态的总数是4。因此，
% 这是一个介于1和4之间的值数组，例如：[1,1,1,1,2,2,2,3,3,3,4,4,4,4,1,1,1]，
% 说明了特征中每个样本的“真”状态标签。
% 状态1=S1声音
% 状态2=收缩
% 状态3=S2声音
% 状态4=舒张
%
% This code was developed by David Springer for comparison purposes in the
% paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
% where a novel segmentation approach is compared to the paper by Schmidt
% et al:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a 
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
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



if(nargin<5)
    figures = false;
end

states = zeros(length(envelope),1);


%% Timing durations from Schmidt:
%%施密特的计时持续时间：
mean_S1 = 0.122*samplingFrequency;
std_S1 = 0.022*samplingFrequency;
mean_S2 = 0.092*samplingFrequency;
std_S2 = 0.022*samplingFrequency;

%% Setting the duration from each R-peak to (R-peak+mean_S1) as the first state:
%%将从每个R-峰值到（R-峰值+平均值）的持续时间设置为第一个状态：
% The R-peak in the ECG coincides with the start of the S1 sound (A. G.
% Tilkian and M. B. Conover, Understanding heart sounds and murmurs: with
% an introduction to lung sounds, 4th ed. Saunders, 2001.)
% Therefore, the duration from each R-peak to the mean_S1 sound duration
% later were labelled as the "true" positions of the S1 sounds:
% 心电图中的R-峰与S1音的开始相吻合（A.G.Tilkian和M.B.Conover，《了解心音和杂音：肺音导论》，第4版，桑德斯，2001年）。
% 因此，从每个R-峰到随后的平均S1音持续时间的持续时间被标记为S1音的“真实”位置：
for i = 1: length(s1_positions)
    %Set an upper bound, incase the window extends over the length of the
    %signal:
    %设置上限，以防窗口超出信号长度：
    upper_bound = round(min(length(states), s1_positions(i) + mean_S1));
    
    %Set the states between the start of the R peak and the upper bound as
    %state 1:
    %将R峰值开始和上限之间的状态设置为状态1：
    states(max([1,s1_positions(i)]):min([upper_bound,length(states)])) = 1;
end

%% Set S2 as state 3 depending on position of end T-wave peak in ECG:
%%根据心电图T波末端峰值的位置，将S2设为状态3：
% The second heart sound occurs at approximately the same time as the
% end-T-wave (A. G. Tilkian and M. B. Conover, Understanding heart sounds
% and murmurs: with an introduction to lung sounds, 4th ed. Saunders, 2001.)
% Therefore, for each end-T-wave, find the peak in the envelope around the
% end-T-wave, setting a window centered on this peak as the second heart
% sound state:

%第二个心音出现的时间与T波末端大致相同（A.G.Tilkian和M.B.Conover，《理解心音和杂音：肺音导论》，第4版，桑德斯，2001年）
%因此，对于每个T波末端，找出T波末端包络线的峰值，将以该峰值为中心的窗口设置为第二个心音状态：

for i = 1: length(s2_positions)
    
    %find search window of envelope:
    %T-end +- mean+1sd
    %Set upper and lower bounds, to avoid errors of searching outside size
    %of the signal
    
    %查找信封的搜索窗口：
    %T端+-平均值+1sd
    %设置上限和下限，以避免搜索超出信号大小的错误
    lower_bound = max([s2_positions(i) - floor((mean_S2 + std_S2)),1]);
    upper_bound = min(length(states), ceil(s2_positions(i) + floor(mean_S2 + std_S2)));
    search_window = envelope(lower_bound:upper_bound).*(states(lower_bound:upper_bound)~=1);
    
    % Find the maximum value of the envelope in the search window:
    %在搜索窗口中查找信封的最大值：
    [~, S2_index] = max(search_window);
    
    %Find the actual index in the envelope of the maximum peak:
    %Make sure this has a max value of the length of the signal:
    %在最大峰值包络中找到实际指数：
    %确保信号长度的最大值为：
    S2_index = min(length(states),lower_bound+ S2_index-1);
    
    %Set the states to state 3, centered on the S2 peak, +- 1/2 of the
    %expected S2 sound duration. Again, making sure it does not try to set a
    %value outside of the length of the signal:
    %将状态设置为状态3，以S2峰值为中心，为预期S2声音持续时间的1/2。同样，确保它不会试图设置超出信号长度的值：
    upper_bound = min(length(states), ceil(S2_index +((mean_S2)/2)));
    states(max([ceil(S2_index - ((mean_S2)/2)),1]):upper_bound) = 3;
    
    %Set the spaces between state 3 and the next R peak as state 4:
    %将状态3和下一个R峰值之间的空间设置为状态4
    if(i<=length(s2_positions))
        %We need to find the next R peak after this S2 sound
        %So, subtract the position of this S2 from the S1 positions
        %我们需要找出这个S2声音后的下一个R峰所以，从S1位置减去S2的位置
        diffs = (s1_positions - s2_positions(i));
        %Exclude those that are negative (meaning before this S2 occured)
        %by setting them to infinity. They are then excluded when finding
        %the minumum later
        %通过将其设置为无穷大，排除那些负数（即在S2发生之前）。当稍后找到最小值时，它们被排除在外
        diffs(diffs<0) = inf;
        
        %If the array is empty, then no S1s after this S2, so set to end of
        %signal:
        %如果阵列为空，则在该S2之后没有S1s，因此设置为信号结束：
        
        if(isempty(diffs<inf))
            end_pos = length(states);
        else
            %else, send the end position to the minimum diff -1
            %否则，将结束位置发送到最小diff-1
            [~, index] = min(diffs);
            end_pos = s1_positions(index) -1;
        end
        states(ceil(S2_index +((mean_S2 +(0*std_S2))/2)):end_pos) = 4;
    end
end




%% Setting the first and last sections of the signal
%%设置信号的第一段和最后一段
% As all states are derived from either R-peak or end-T-wave locations, the first affirmed state
% in the signal will always be state 1 or state 3. Therefore, until this state, the
% first state should always be set to 4 or 2:
%由于所有的状态都是从R-峰值或T-端波位置导出的，信号中的第一个确认状态总是状态1或状态3。
%因此，在此状态之前，第一个状态应始终设置为4或2：

%Find the first step up:
first_location_of_definite_state = find(states ~= 0, 1)-1;

if(first_location_of_definite_state > 1)
    
    if(states(first_location_of_definite_state + 1) == 1)
        states(1:first_location_of_definite_state) = 4;
    end
    
    if(states(first_location_of_definite_state + 1) == 3)
        states(1:first_location_of_definite_state) = 2;
    end
    
end


% Find the last step down:
last_location_of_definite_state = find(states ~= 0, 1,'last');

if(last_location_of_definite_state > 1)
    
    if(states(last_location_of_definite_state) == 1)
        states(last_location_of_definite_state:end) = 2;
    end
    
    if(states(last_location_of_definite_state) == 3)
        states(last_location_of_definite_state:end) = 4;
    end
    
end


states(length(envelope)+1 : end) = [];


%Set everywhere else as state 2:
states(states == 0) = 2;


%% Plotting figures
if(figures)
    figure('Name','Envelope and labelled states');
    plot(envelope);
    hold on;
    plot(states,'r');
    legend('Envelope', 'States');
    pause();
end



