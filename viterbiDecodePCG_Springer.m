% function [delta, psi, qt] = viterbiDecodePCG_Springer(observation_sequence, pi_vector, b_matrix, total_obs_distribution, heartrate, systolic_time, Fs, figures)
%
% This function calculates the delta, psi and qt matrices associated with
% the Viterbi decoding algorithm from:
%此函数根据以下公式计算与维特比解码算法相关的增量、psi和qt矩阵：
% L. R. Rabiner, "A tutorial on hidden Markov models and selected
% applications in speech recognition," Proc. IEEE, vol. 77, no. 2, pp.0
% 257-286, Feb. 1989.
% using equations 32a - 35, and equations 68 - 69 to include duration
% dependancy of the states.
%
% This decoding is performed after the observation probabilities have been
% derived from the logistic regression model of Springer et al:
%在从Springer等人的logistic回归模型推导出观察概率后，进行解码：
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
% Further, this function is extended to allow the duration distributions to extend
% past the beginning and end of the sequence. Without this, the label
% sequence has to start and stop with an "entire" state duration being
% fulfilled. This extension takes away that requirement, by allowing the
% duration distributions to extend past the beginning and end, but only
% considering the observations within the sequence for emission probability
% estimation. More detail can be found in the publication by Springer et
% al., mentioned above.
%此外，该函数被扩展以允许持续时间分布扩展到序列的开始和结束。
%否则，标签序列必须以“整个”状态持续时间开始和停止。
%这种扩展消除了这一要求，允许持续时间分布扩展到超过开始和结束，
%但只考虑排放概率估计序列内的观测值。更多细节可以在上面提到的Springer等人的出版物中找到。
%% Inputs:
% observation_sequence: The observed features 
% pi_vector: the array of initial state probabilities, dervived from
% "trainSpringerSegmentationAlgorithm".
% b_matrix: the observation probabilities, dervived from
% "trainSpringerSegmentationAlgorithm".
% heartrate: the heart rate of the PCG, extracted using
% "getHeartRateSchmidt"
% systolic_time: the duration of systole, extracted using
% "getHeartRateSchmidt"
% Fs: the sampling frequency of the observation_sequence
% figures: optional boolean variable to show figures
%
% 观测序列：观察到的特征
% pi向量：由“trainSpringerSegmentationAlgorithm”导出的初始状态概率数组。
% b矩阵：观测概率，来源于“trainSpringerSegmentationAlgorithm”。
% 心率：PCG的心率，使用“getHeartRateSchmidt”提取
% 收缩时间：收缩持续时间，使用“GetHeartrateSmidt”提取
% Fs：观察序列图形的采样频率：可选布尔变量显示图形
%
%% Outputs:
% logistic_regression_B_matrix:
% pi_vector:
% total_obs_distribution:
% As Springer et al's algorithm is a duration dependant HMM, there is no
% need to calculate the A_matrix, as the transition between states is only
% dependant on the state durations.
%
% 逻辑回归B矩阵
% pi向量
% total_obs_distribution:
% 由于Springer等人的算法是一个依赖于持续时间的HMM，因此不需要计算a_矩阵，因为状态之间的转换只依赖于状态持续时间。
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

function [delta, psi, qt] = viterbiDecodePCG_Springer(observation_sequence, pi_vector, b_matrix, total_obs_distribution, heartrate, systolic_time, Fs,figures)

if nargin < 8
    figures = false;
end

%% Preliminary 初始步骤
springer_options = default_Springer_HSMM_options;

T = length(observation_sequence);
N = 4; % Number of states

% Setting the maximum duration of a single state. This is set to an entire
% heart cycle:
% 设置单个状态的最大持续时间。这将设置为整个心脏循环：
max_duration_D = round((1*(60/heartrate))*Fs);

% Initialising the variables that are needed to find the optimal state path along
% the observation sequence.
% 初始化沿观测序列寻找最优状态路径所需的变量。
% delta_t(j), as defined on page 264 of Rabiner, is the best score (highest
% probability) along a single path, at time t, which accounts for the first
% t observations and ends in State s_j. In this case, the length of the
% matrix is extended by max_duration_D samples, in order to allow the use
% of the extended Viterbi algortithm:
% delta_t（j）如Rabiner第264页所定义，是在时间t沿单个路径的最佳分数（最高概率），
% 它说明了第一个t观测值并以状态s_j结束。在这种情况下，矩阵的长度通过最大持续时间D样本扩展，
% 以便允许使用扩展的Viterbi算法：
delta = ones(T+ max_duration_D-1,N)*-inf;

% The argument that maximises the transition between states (this is
% basically the previous state that had the highest transition probability
% to the current state) is tracked using the psi variable.
% 使用psi变量跟踪最大化状态之间转换的论点（这基本上是以前的状态，它具有最高的过渡概率到当前状态）。
psi = zeros(T+ max_duration_D-1,N);

% An additional variable, that is not included on page 264 or Rabiner, is
% the state duration that maximises the delta variable. This is essential
% for the duration dependant HMM.
% 第264页或Rabiner中未包含的另一个变量是使delta变量最大化的状态持续时间。这对于依赖于持续时间的HMM来说是必不可少的。
psi_duration =zeros(T + max_duration_D-1,N);

%% Setting up observation probs 设置观察问题
observation_probs = zeros(T,N);

for n = 1:N
    
    % MLR gives P(state|obs)
    % Therefore, need Bayes to get P(o|state)
    % P(o|state) = P(state|obs) * P(obs) / P(states)
    % Where p(obs) is derived from a MVN distribution from all
    % obserbations, and p(states) is taken from the pi_vector:
    pihat = mnrval(cell2mat(b_matrix(n)),observation_sequence(:,:));
    
    for t = 1:T
        
        Po_correction = mvnpdf(observation_sequence(t,:),cell2mat(total_obs_distribution(1)),cell2mat(total_obs_distribution(2)));
        
        %When saving the coefficients from the logistic
        %regression, it orders them P(class 1) then P(class 2). When
        %training, I label the classes as 0 and 1, so the
        %correct probability would be pihat(2).
        % 当保存logistic回归的系数时，它将它们排序为P（类1），然后是P（类2）。
        % 训练时，我将类标记为0和1，因此正确的概率是pihat（2）。
        
        observation_probs(t,n) = (pihat(t,2)*Po_correction)/pi_vector(n);
        
    end
end

%% Setting up state duration probabilities, using Gaussian distributions:
%使用高斯分布设置状态持续时间概率：
[d_distributions, max_S1, min_S1, max_S2, min_S2, max_systole, min_systole, max_diastole, min_diastole] = get_duration_distributions(heartrate,systolic_time);



duration_probs = zeros(N,3*Fs);
duration_sum = zeros(N,1);
for state_j = 1:N
    for d = 1:max_duration_D
        if(state_j == 1)
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_S1 || d > max_S1)
                duration_probs(state_j,d)= realmin;
            end
            
            
        elseif(state_j==3)
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_S2 || d > max_S2)
                duration_probs(state_j,d)= realmin;
            end
            
            
        elseif(state_j==2)
            
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_systole|| d > max_systole)
                duration_probs(state_j,d)= realmin;
            end
            
            
        elseif (state_j==4)
            
            duration_probs(state_j,d) = mvnpdf(d,cell2mat(d_distributions(state_j,1)),cell2mat(d_distributions(state_j,2)));
            
            if(d < min_diastole ||d > max_diastole)
                duration_probs(state_j,d)= realmin;
            end
        end
    end
    duration_sum(state_j) = sum(duration_probs(state_j,:));
end


if(length(duration_probs)>3*Fs)
    duration_probs(:,(3*Fs+1):end) = [];
end

if(figures)
    figure('Name', 'Duration probabilities');
    plot(duration_probs(1,:)./ duration_sum(1),'Linewidth',2);
    hold on;
    plot(duration_probs(2,:)./ duration_sum(2),'r','Linewidth',2);
    hold on;
    plot(duration_probs(3,:)./ duration_sum(3),'g','Linewidth',2);
    hold on;
    plot(duration_probs(4,:)./ duration_sum(4),'k','Linewidth',2);
    hold on;
    legend('S1 Duration','Systolic Duration','S2 Duration','Diastolic Duration');
    pause();
end
%% Perform the actual Viterbi Recursion:
%执行实际的Viterbi递归：


qt = zeros(1,length(delta));
%% Initialisation Step
% 初始化步骤

%Equation 32a and 69a, but leave out the probability of being in
%state i for only 1 sample, as the state could have started before time t =
%0.

delta(1,:) = log(pi_vector) + log(observation_probs(1,:)); %first value is the probability of intially being in each state * probability of observation 1 coming from each state

%Equation 32b
psi(1,:) = -1;


% The state duration probabilities are now used.
%Change the a_matrix to have zeros along the diagonal, therefore, only
%relying on the duration probabilities and observation probabilities to
%influence change in states:
%This would only be valid in sequences where the transition between states
%follows a distinct order.
% 现在使用状态持续时间概率。因此，将a_矩阵更改为沿对角线有零，仅依赖持续时间概率
% 和观察概率来影响状态变化：这仅在状态之间的转换遵循不同顺序的序列中有效。
a_matrix = [0,1,0,0;0 0 1 0; 0 0 0 1;1 0 0 0];


%% Run the core Viterbi algorith 运行核心维特比算法

if(springer_options.use_mex)
    
    %% Run Mex code 运行Mex代码
    % Ensure you have run the mex viterbi_PhysChallenge.c code on the
    % native machine before running this:
    %在运行此操作之前，请确保已在本机计算机上运行mex viterbi\u PhysChallenge.c代码：
    [delta, psi, psi_duration] = viterbi_Springer(N,T,a_matrix,max_duration_D,delta,observation_probs,duration_probs,psi, duration_sum);
    
    
else
    
    %% Recursion 递归
    
    %% The Extended Viterbi algorithm: 扩展维特比算法：
    
    %Equations 33a and 33b and 69a, b, c etc:
    %again, ommitting the p(d), as state could have started before t = 1
    
    % This implementation extends the standard implementation of the
    % duration-dependant Viterbi algorithm by allowing the durations to
    % extend beyond the start and end of the time series, thereby allowing
    % states to "start" and "stop" outside of the recorded signal. This
    % addresses the issue of partial states at the beginning and end of the
    % signal being labelled as the incorrect state. For instance, a
    % short-duration diastole at the beginning of a signal looks a lot like
    % systole, and can lead to labelling errors.
    %该操作通过允许持续时间扩展到时间序列的开始和结束之外，
    %从而允许状态在所记录的信号之外“开始”和“停止”，
    %从而扩展了持续时间相关的Viterbi算法的标准实现。
    %这解决了信号开头和结尾部分状态被标记为错误状态的问题。
    %例如，在信号开始时的短暂舒张期看起来很像收缩期，可能导致标记错误。
    
    % t spans input 2 to T + max_duration_D:
    
    
    for t = 2:T+ max_duration_D-1
        for j = 1:N
            for d = 1:1:max_duration_D
                
                
                %The start of the analysis window, which is the current time
                %step, minus d (the time horizon we are currently looking back),
                %plus 1. The analysis window can be seen to be starting one
                %step back each time the variable d is increased.
                % This is clamped to 1 if extending past the start of the
                % record, and T-1 is extending past the end of the record:
                start_t = t - d;
                if(start_t<1)
                    start_t = 1;
                end
                if(start_t > T-1)
                    start_t = T-1;
                end
                
                %The end of the analysis window, which is the current time
                %step, unless the time has gone past T, the end of the record, in
                %which case it is truncated to T. This allows the analysis
                %window to extend past the end of the record, so that the
                %timing durations of the states do not have to "end" at the end
                %of the record.
                end_t = t;
                if(t>T)
                    end_t = T;
                end
                
                
                %Find the max_delta and index of that from the previous step
                %and the transition to the current step:
                %This is the first half of the expression of equation 33a from
                %Rabiner:
                [max_delta, max_index] = max(delta(start_t,:)+log(a_matrix(:,j))');
                               
                
                %Find the normalised probabilities of the observations over the
                %analysis window:
                probs = prod(observation_probs(start_t:end_t,j));
                
                
                %Find the normalised probabilities of the observations at only
                %the time point at the start of the time window:
                
                if(probs ==0)
                    probs = realmin;
                end
                emission_probs = log(probs);
                
                
                %Keep a running total of the emmission probabilities as the
                %start point of the time window is moved back one step at a
                %time. This is the probability of seeing all the observations
                %in the analysis window in state j:
                
                if(emission_probs == 0 || isnan(emission_probs))
                    emission_probs =realmin;
                end
                
                
                %Find the total probability of transitioning from the last
                %state to this one, with the observations and being in the same
                %state for the analysis window. This is the duration-dependant
                %variation of equation 33a from Rabiner:
                %                 fprintf('log((duration_probs(j,d)./duration_sum(j))):%d\n',log((duration_probs(j,d)./duration_sum(j))));
                delta_temp = max_delta + (emission_probs)+ log((duration_probs(j,d)./duration_sum(j)));
                
                
                %Unlike equation 33a from Rabiner, the maximum delta could come
                %from multiple d values, or from multiple size of the analysis
                %window. Therefore, only keep the maximum delta value over the
                %entire analysis window:
                %If this probability is greater than the last greatest,
                %update the delta matrix and the time duration variable:
                
                
                if(delta_temp>delta(t,j))
                    delta(t,j) = delta_temp;
                    psi(t,j) = max_index;
                    psi_duration(t,j) = d;
                end
                
            end
        end
    end
end


%% Termination 终止

% For the extended case, need to find max prob after end of actual
% sequence:
% 对于扩展情况，需要在实际序列结束后找到max prob：

% Find just the delta after the end of the actual signal
% 找到实际信号结束后的增量
temp_delta = delta(T+1:end,:);
% Find the maximum value in this section, and which state it is in:
% 找出本节中的最大值，以及它所处的状态：
[~, idx] = max(temp_delta(:));
[pos, ~] = ind2sub(size(temp_delta), idx);

% Change this position to the real position in delta matrix:
% 将此位置更改为增量矩阵中的实际位置：
pos = pos+T;

%1) Find the last most probable state
%2) From the psi matrix, find the most likely preceding state
%3) Find the duration of the last state from the psi_duration matrix
%4) From the onset to the offset of this state, set to the most likely state
%5) Repeat steps 2 - 5 until reached the beginning of the signal
% 找到最后一个最可能的状态
% 从psi矩阵中，找出最可能的前一状态 
% 从psi\U持续时间矩阵中查找最后一个状态的持续时间
% 从该状态的开始到偏移，设置为最可能的状态
% 重复步骤2-5，直到到达信号的开头


%The initial steps 1-4 are equation 34b in Rabiner. 1) finds P*, the most
%likely last state in the sequence, 2) finds the state that precedes the
%last most likely state, 3) finds the onset in time of the last state
%(included due to the duration-dependancy) and 4) sets the most likely last
%state to the q_t variable.
% 初始步骤1-4是Rabiner中的等式34b。1） 查找P*，序列中最可能的最后一个状态，
% 2）查找最可能出现的状态，3）查找最后一个状态的开始时间（由于持续时间依赖性而包含），
% 4）将最可能的最后一个状态设置为q_ut变量。

%1)
[~, state] = max(delta(pos,:),[],2);

%2)
offset = pos;
preceding_state = psi(offset,state);

%3)
% state_duration = psi_duration(offset, state);
onset = offset - psi_duration(offset,state)+1;

%4)
qt(onset:offset) = state;

%The state is then updated to the preceding state, found above, which must
%end when the last most likely state started in the observation sequence:
% 然后将状态更新为前面的状态（见上文），该状态必须在观察序列中最后一个最可能的状态开始时结束：
state = preceding_state;

count = 0;
%While the onset of the state is larger than the maximum duration
%specified:
% 当状态开始时间大于指定的最大持续时间时：
while(onset > 2)
    
    %2)
    offset = onset-1;
    %     offset_array(offset,1) = inf;
    preceding_state = psi(offset,state);
    %     offset_array(offset,2) = preceding_state;
    
    
    %3)
    %     state_duration = psi_duration(offset, state);
    onset = offset - psi_duration(offset,state)+1;
    
    %4)
    %     offset_array(onset:offset,3) = state;
    
    if(onset<2)
        onset = 1;
    end
    qt(onset:offset) = state;
    state = preceding_state;
    count = count +1;
    
    if(count> 1000)
        break;
    end
end

qt = qt(1:T);


