% function [B_matrix, pi_vector, total_obs_distribution] = trainBandPiMatricesSpringer(state_observation_values)
%
% Train the B matrix and pi vector for the Springer HMM.ΪSpringer HMMѵ��B�����pi������
% The pi vector is the initial state probability, while the B matrix are
% the observation probabilities. In the case of Springer's algorith, the
% observation probabilities are based on a logistic regression-based
% probabilities. 
%pi�����ǳ�ʼ״̬���ʣ���B�����ǹ۲���ʡ���Springer�㷨������£��۲���ʻ��ڻ���logistic�ع�ĸ��ʡ�
%% ʹ���߼��ع�ó�B����ѵ��HSMM
%% Inputs:
% state_observation_values: an Nx4 cell array of observation values from
% each of N PCG signals for each (of 4) state. Within each cell is a KxJ
% double array, where K is the number of samples from that state in the PCG
% and J is the number of feature vectors extracted from the PCG.
%״̬�۲�ֵ��ÿ������4����״̬��N��PCG�źŵĹ۲�ֵ��Nx4��Ԫ���С�ÿ����Ԫ����һ��KxJ˫���飬
%����K��PCG�и�״̬����������J�Ǵ�PCG����ȡ������������������

%% Outputs:
% The B_matrix and pi arrays for an HMM - as Springer et al's algorithm is a
% duration dependant HMM, there is no need to calculate the A_matrix, as
% the transition between states is only dependant on the state durations.
% total_obs_distribution:
%������ɷ�ģ�͵�B�����pi����-����Springer���˵��㷨��һ�������ʱ����ص�HMM��
%��˲���Ҫ����a_������Ϊ״̬֮���ת��ֻ������״̬����ʱ�䡣�ֲܷ���

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
%pi��������ֵ������ʼ״̬���ʣ�ȡ����ÿ��PCG�����ʺ�ÿ�����ߵ���������ʱ�䡣
% ��ˣ��������û���������pi_����������ͨ��������״̬����Ϊ��ͬ�ĸ������򻯣�


pi_vector = [0.25,0.25,0.25,0.25];

%% Train the logistic regression-based B_matrix:


% Initialise the B_matrix as a 1x4 cell array. This is to hold the
% coefficients of the trained logisitic regression model for each state.
% ��B_�����ʼ��Ϊ1x4��Ԫ���С�����Ϊ�˱���ÿ��״̬��ѵ�����߼��ع�ģ�͵�ϵ����
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

% Ϊ��ʹ�ñ�Ҷ˹��ʽ��logistic�ع鵼���ĸ��ʣ�������Ҫ�õ�������ѵ�����ݼ��п���
% �ض��۲�ֵ�ĸ��ʡ����ǡ�����۲����С���ÿ��״̬��ƽ��ֵ��Э�������£�

total_observation_sequence = vertcat(statei_values{1}, statei_values{2}, statei_values{3}, statei_values{4});
total_obs_distribution = cell(2,1);
total_obs_distribution{1} = mean(total_observation_sequence);
total_obs_distribution{2} = cov(total_observation_sequence);


for state = 1: number_of_states
    
    % Randomly select indices of samples from the other states not being 
    % learnt, in order to balance the two data sets. The code below ensures
    % that if class 1 is being learnt vs the rest, the number of the rest =
    % the number of class 1, evenly split across all other classes
    % ���ѡȡ����δѧϰ״̬������ָ�꣬��ƽ���������ݡ�����Ĵ���ȷ���ˣ�
    % �����1�����������������ѧϰ����ʣ�������=��1���������ƽ�����������������
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
    %����ζ������һ��״̬û���㹻��������ƥ������ѵ������ĳ��ȣ�
    
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
            %ȷ��ֻΪ��״̬ѡ��n-1��*3*length_per_����_state samples��Ϊ��״̬����ȷ��������ƽ��ģ�
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

