load('PCGdata');
  
% ���������
B=randperm(length(PCGdata(1,:)));
data=PCGdata(:,B(1:length(PCGdata(1,:)))); 
% ����ѵ�����Ͳ��Լ�
A = (1:515);
C = (516:644);
traindata=data(:,A(1:515)); 
testdata = data(:,C(1:129)); 
    
train = traindata{1,1};

[Mdl,FitInfo] = fitclinear(X,Ystats);