function [ ss ] = feature_extraction(n)
%addpath(genpath(' D:\MATLAB\R2017a\toolbox\voicebox'))
load('new_example_data');
[x,]=new_data{1,n};
% 求梅尔倒谱系数，并压缩所有帧数为1X16的行向量
ccc = mfcc_m(x);
cc = std(ccc,0,1);
% 求过零率，一个数值
count = zero_crossings(x);
% 求能量均方根，一个数值
RMSE = rms(x);
% 求经验模态分解，并对矩阵求均值，将其压缩为一个数值
y = emd(x,'fix',10);
yy = mean(mean(y));
%构建特征向量
ss = [count,RMSE,cc,yy];
 end