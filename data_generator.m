%%
% 数据预处理，包括剔除无效数据，随机化，划分训练集，测试集
clear
clc
% 导入数据，去除无效数据
load('PCGdata')
i = 1;
for j = 1:792
if isempty(PCGdata{1,j})
        a(i)=j;
        i=i+1;
end
end

%PCGdata(:,a)=[];

% 随机化处理
b=randperm(length(PCGdata(1,:)));
data=PCGdata(:,b(1:length(PCGdata(1,:)))); 
% 分配训练集和测试集
a = (1:515);
c = (516:644);
traindata=data(:,a(1:515)); 
testdata = data(:,c(1:129)); 



