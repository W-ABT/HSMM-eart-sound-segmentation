%%
% ����Ԥ���������޳���Ч���ݣ������������ѵ���������Լ�
clear
clc
% �������ݣ�ȥ����Ч����
load('PCGdata')
i = 1;
for j = 1:792
if isempty(PCGdata{1,j})
        a(i)=j;
        i=i+1;
end
end

%PCGdata(:,a)=[];

% ���������
b=randperm(length(PCGdata(1,:)));
data=PCGdata(:,b(1:length(PCGdata(1,:)))); 
% ����ѵ�����Ͳ��Լ�
a = (1:515);
c = (516:644);
traindata=data(:,a(1:515)); 
testdata = data(:,c(1:129)); 



