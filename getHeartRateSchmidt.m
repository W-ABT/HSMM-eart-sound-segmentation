% function [heartRate systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)
%
% Derive the heart rate and the sytolic time interval from a PCG recording.
% This is used in the duration-dependant HMM-based segmentation of the PCG
% recording.
%从PCG记录中得出心率和综合时间间隔。这用于基于持续时间的基于HMM的PCG记录分割
%
% This method is based on analysis of the autocorrelation function, and the
% positions of the peaks therein.
%该方法基于对自相关函数及其峰值位置的分析。
%
% This code is derived from the paper:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a 
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% INPUTS:
% audio_data: The raw audio data from the PCG recording  音频数据：来自PCG记录的原始音频数据
% Fs: the sampling frequency of the audio recording      Fs：录音的采样频率
% figures: optional boolean to display figures           图片：显示图形的可选布尔值
%
%% OUTPUTS:
% heartRate: the heart rate of the PCG in beats per minute            心率：PCG的心率，单位为每分钟心跳数
% systolicTimeInterval: the duration of systole, as derived from the
% autocorrelation function, in seconds                  收缩间期：由自相关函数导出的收缩持续时间，以秒为单位
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

function [heartRate, systolicTimeInterval] = getHeartRateSchmidt(audio_data, Fs, figures)

if nargin < 3
    figures = false;
end



%输入参数必须为三类，因为输入参数必须包括audio_data（音频数据），Fs（采样频率），figures（图片）
%% Get heatrate:
% From Schmidt:
% "The duration of the heart cycle is estimated as the time from lag zero
% to the highest peaks between 500 and 2000 ms in the resulting
% autocorrelation"
% This is performed after filtering and spike removal:
%“心脏周期的持续时间估计为自相关中从滞后零点到500到2000毫秒之间的最高峰值的时间”，这是在过滤和去除尖峰之后执行的：

%% 25-400Hz 4th order Butterworth band pass                   25-400Hz四阶巴特沃斯带通
audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs, false);
audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs);

%% Spike removal from the original paper:                     去除原文件中的噪音
audio_data = schmidt_spike_removal(audio_data,Fs);

%% Find the homomorphic envelope                              找到同态包络
homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs);

%% Find the autocorrelation:                                  求自相关：
y=homomorphic_envelope-mean(homomorphic_envelope);
[c] = xcorr(y,'coeff');
signal_autocorrelation = c(length(homomorphic_envelope)+1:end);

min_index = 0.5*Fs;
max_index = 2*Fs;

[~, index] = max(signal_autocorrelation(min_index:max_index));
true_index = index+min_index-1;

heartRate = 60/(true_index/Fs);


%% Find the systolic time interval:                          找出收缩时间间隔
% From Schmidt: "The systolic duration is defined as the time from lag zero
% to the highest peak in the interval between 200 ms and half of the heart
% cycle duration"
%收缩持续时间是指从滞后零点到200毫秒到心脏周期一半时间间隔内的最高峰值的时间。

max_sys_duration = round(((60/heartRate)*Fs)/2);
min_sys_duration = round(0.2*Fs);

[~, pos] = max(signal_autocorrelation(min_sys_duration:max_sys_duration));
systolicTimeInterval = (min_sys_duration+pos)/Fs;


if(figures)
    figure('Name', 'Heart rate calculation figure');
    plot(signal_autocorrelation);
    hold on;
    plot(true_index, signal_autocorrelation(true_index),'ro');
    plot((min_sys_duration+pos), signal_autocorrelation((min_sys_duration+pos)), 'mo');
    xlabel('Samples');
    legend('Autocorrelation', 'Position of max peak used to calculate HR', 'Position of max peak within systolic interval');
end

