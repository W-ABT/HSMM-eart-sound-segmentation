% function [B_matrix, pi_vector, total_obs_distribution] = trainBandPiMatricesSpringer(state_observation_values)
%
% Train the B matrix and pi vector for the Springer HMM.为Springer HMM训练B矩阵和pi向量。
% The pi vector is the initial state probability, while the B matrix are
% the observation probabilities. In the case of Springer's algorith, the
% observation probabilities are based on a logistic regression-based
% probabilities. 
%pi向量是初始状态概率，而B矩阵是观测概率。在Springer算法的情况下，观察概率基于基于logistic回归的概率。
%% 使用逻辑回归得出B矩阵，训练HSMM
%% Inputs:
% state_observation_values: an Nx4 cell array of observation values from
% each of N PCG signals for each (of 4) state. Within each cell is a KxJ
% double array, where K is the number of samples from that state in the PCG
% and J is the number of feature vectors extracted from the PCG.
%状态观察值：每个（共4个）状态的N个PCG信号的观察值的Nx4单元阵列。每个单元内有一个KxJ双数组，
%其中K是PCG中该状态的样本数，J是从PCG中提取的特征向量的数量。

%% Outputs:
% The B_matrix and pi arrays for an HMM - as Springer et al's algorithm is a
% duration dependant HMM, there is no need to calculate the A_matrix, as
% the transition between states is only dependant on the state durations.
% total_obs_distribution:
%隐马尔可夫模型的B矩阵和pi数组-由于Springer等人的算法是一个与持续时间相关的HMM，
%因此不需要计算a_矩阵，因为状态之间的转换只依赖于状态持续时间。总分布：

% Developed by David Springer for the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
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

function [B_matrix, pi_vector, total_obs_distribution] = trainBandPiMatricesSpringer(state_observation_values)

%% Prelim

number_of_states = 4;

%% Set pi_vector
% The true value of the pi vector, which are the initial state
% probabilities, are dependant on the heart rate of each PCG, and the
% individual sound duration for each patient. Therefore, instead of setting
% a patient-dependant pi_vector, simplify by setting all states as equally
% probable:
%pi向量的真值，即初始状态概率，取决于每个PCG的心率和每个患者的声音持续时间。
% 因此，与其设置患者依赖的pi_向量，不如通过将所有状态设置为相同的概率来简化：


pi_vector = [0.25,0.25,0.25,0.25];

%% Train the logistic regression-based B_matrix:


% Initialise the B_matrix as a 1x4 cell array. This is to hold the
% coefficients of the trained logisitic regression model for each state.
% 将B_矩阵初始化为1x4单元阵列。这是为了保存每个状态下训练的逻辑回归模型的系数。
B_matrix = cell(1,number_of_states);

statei_values = cell(number_of_states,1);

for PCGi = 1: length(state_observation_values)
        
    statei_values{1} = vertcat(statei_values{1},state_observation_values{PCGi,1});
    statei_values{2} = vertcat(statei_values{2},state_observation_values{PCGi,2});
    statei_values{3} = vertcat(statei_values{3},state_observation_values{PCGi,3});
    statei_values{4} = vertcat(statei_values{4},state_observation_values{PCGi,4});
    
end


% In order to use Bayes' formula with the logistic regression derived
% probabilities, we need to get the probability of seeing a specific
% observation in the total training data set. This is the
% 'total_observation_sequence', and the mean and covariance for each state
% is found:

% 为了使用贝叶斯公式和logistic回归导出的概率，我们需要得到在整个训练数据集中看到
% 特定观测值的概率。这是“总体观察序列”，每个状态的平均值和协方差如下：

total_observation_sequence = vertcat(statei_values{1}, statei_values{2}, statei_values{3}, statei_values{4});
total_obs_distribution = cell(2,1);
total_obs_distribution{1} = mean(total_observation_sequence);
total_obs_distribution{2} = cov(total_observation_sequence);


for state = 1: number_of_states
    
    % Randomly select indices of samples from the other states not being 
    % learnt, in order to balance the two data sets. The code below ensures
    % that if class 1 is being learnt vs the rest, the number of the rest =
    % the number of class 1, evenly split across all other classes
    % 随机选取其他未学习状态样本的指标，以平衡两组数据。下面的代码确保了，
    % 如果第1类与其他类相比正在学习，则剩余的数量=第1类的数量，平均分配给所有其他类
    length_of_state_samples = length(statei_values{state});
    
    % Number of samples required from each of the other states:
    length_per_other_state = floor(length_of_state_samples/(number_of_states-1));
    
    
    %If the length of the main class / (num states - 1) >
    %length(shortest other class), then only select
    %length(shortect other class) from the other states,
    %and (3* length) for main class
    min_length_other_class = inf;
    
    for other_state = 1: number_of_states
        samples_in_other_state = length(statei_values{other_state});
        
        if(other_state == state)
        else
            min_length_other_class = min([min_length_other_class, samples_in_other_state]);
        end
    end
    
    %This means there aren't enough samples in one of the
    %states to match the length of the main class being
    %trained:
    %这意味着其中一个状态没有足够的样本来匹配所培训的主类的长度：
    
    if( length_per_other_state > min_length_other_class)
        length_per_other_state = min_length_other_class;
    end
    
    training_data = cell(2,1);
    
    for other_state = 1: number_of_states
        samples_in_other_state = length(statei_values{other_state});
                
        if(other_state == state)
            %Make sure you only choose (n-1)*3 *
            %length_per_other_state samples for the main
            %state, to ensure that the sets are balanced:
            %确保只为主状态选择（n-1）*3*length_per_其他_state samples作为主状态，以确保集合是平衡的：
            indices = randperm(samples_in_other_state,length_per_other_state*(number_of_states-1));
            training_data{1} = statei_values{other_state}(indices,:);
        else
                       
            indices = randperm(samples_in_other_state,length_per_other_state);
            state_data = statei_values{other_state}(indices,:);
            training_data{2} = vertcat(training_data{2}, state_data);
            
        end
    end
    
    % Label all the data:
    labels = ones(length(training_data{1}) + length(training_data{2}),1);
    labels(1:length(training_data{1})) = 2;
    
    % Train the logisitic regression model for this state:
    all_data = [training_data{1};training_data{2}];
    [B,~,~] = mnrfit(all_data,labels);
    B_matrix{state} = B;
end

