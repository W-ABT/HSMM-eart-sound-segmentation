%% 
% 脚本，用于生成数据
clear;
clc;
%%
% 默认选项
springer_options = default_Springer_HSMM_options;
load('example_data.mat');
training_indices = [1, 47, 361, 402, 572];
train_recordings = example_data.example_audio_data(training_indices);
train_annotations = example_data.example_annotations(training_indices,:);
[B_matrix, pi_vector, total_obs_distribution] = trainSpringerSegmentationAlgorithm(train_recordings,train_annotations,springer_options.audio_Fs, false);
%%
% 生成数据
test_index = [1:792];

test_recordings = example_data.example_audio_data(test_index);
test_annotations = example_data.example_annotations(test_index,:);
test_binary = example_data.binary_diagnosis(test_index);
numPCGs = length(test_recordings);

for PCGi = 1:numPCGs
    if length(test_recordings{PCGi})>2000
    [assigned_states,PCG_states,PCG_Features] = runSpringerSegmentationAlgorithm(test_recordings{PCGi}, springer_options.audio_Fs, B_matrix, pi_vector, total_obs_distribution, test_annotations,true);
    
    data{1,PCGi} = [test_recordings{PCGi},assigned_states];
    
    for i = 1:length(assigned_states)
        for j =1:4
            if data{1,PCGi}(i,2)==j
                PCGdata{1,PCGi}(i,j)=data{1,PCGi}(i,1);
            end
        end
    end
    PCGdata{2,PCGi}=test_binary{PCGi};
    PCGdata{3,PCGi}=PCGi;
    end
end
