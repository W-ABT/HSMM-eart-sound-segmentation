function [ ss ] = feature_extraction(n)
%addpath(genpath(' D:\MATLAB\R2017a\toolbox\voicebox'))
load('new_example_data');
[x,]=new_data{1,n};
% ��÷������ϵ������ѹ������֡��Ϊ1X16��������
ccc = mfcc_m(x);
cc = std(ccc,0,1);
% ������ʣ�һ����ֵ
count = zero_crossings(x);
% ��������������һ����ֵ
RMSE = rms(x);
% ����ģ̬�ֽ⣬���Ծ������ֵ������ѹ��Ϊһ����ֵ
y = emd(x,'fix',10);
yy = mean(mean(y));
%������������
ss = [count,RMSE,cc,yy];
 end