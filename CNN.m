%%
% CNN算法，训练CNN模型并用来分类
clear
clc
%%
% 导入数据
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')

load('PCGdata');
% 随机化处理
b=randperm(length(PCGdata(1,:)));
data=PCGdata(:,b(1:length(PCGdata(1,:)))); 
% 分配训练集和测试集
a = (1:515);
c = (516:644);
traindata=data(:,a(1:515)); 
testdata = data(:,c(1:129)); 
for i = 1:515
traindata{2,i}=categorical(traindata{2,i});
end
for i = 1:129
testdata{2,i}=categorical(testdata{2,i});
end
train = traindata{1,1};
%%
% 用来设定修改参数的语句如下：
options = trainingOptions('sgdm','MaxEpochs',20,...
    'InitialLearnRate',0.0001);
     
%%
% 在matlab中用来建立网络的语句如下：
layers = [imageInputLayer([length(traindata{1,1}(:,1)) 4 1]);
          convolution2dLayer(5,20);
          reluLayer();
          fullyConnectedLayer(1);
          softmaxLayer();
          classificationLayer()];

%%
% 训练CNN模型
ile = gpuArray(0.0001);
net = trainNetwork(train,traindata{2,1},layers,options);
%%
% 测试集效果

% YTest = classify(convnet,testDigitData);
% TTest = testDigitData.Labels;
% accuracy = sum(YTest == TTest)/numel(TTest);
% disp(accuracy);
