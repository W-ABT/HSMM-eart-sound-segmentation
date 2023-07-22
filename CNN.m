%%
% CNN�㷨��ѵ��CNNģ�Ͳ���������
clear
clc
%%
% ��������
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')

load('PCGdata');
% ���������
b=randperm(length(PCGdata(1,:)));
data=PCGdata(:,b(1:length(PCGdata(1,:)))); 
% ����ѵ�����Ͳ��Լ�
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
% �����趨�޸Ĳ�����������£�
options = trainingOptions('sgdm','MaxEpochs',20,...
    'InitialLearnRate',0.0001);
     
%%
% ��matlab���������������������£�
layers = [imageInputLayer([length(traindata{1,1}(:,1)) 4 1]);
          convolution2dLayer(5,20);
          reluLayer();
          fullyConnectedLayer(1);
          softmaxLayer();
          classificationLayer()];

%%
% ѵ��CNNģ��
ile = gpuArray(0.0001);
net = trainNetwork(train,traindata{2,1},layers,options);
%%
% ���Լ�Ч��

% YTest = classify(convnet,testDigitData);
% TTest = testDigitData.Labels;
% accuracy = sum(YTest == TTest)/numel(TTest);
% disp(accuracy);
