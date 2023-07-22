% function [PCG_Features, featuresFs] = getSpringerPCGFeatures(audio_data, Fs, figures)
%
% Get the features used in the Springer segmentation algorithm. These 
% features include:
% -The homomorphic envelope (as performed in Schmidt et al's paper)
% -The Hilbert envelope
% -A wavelet-based feature
% -A PSD-based feature
% 得到Springer分割算法中用到的特征。这些功能包括：
%-同态包络（如Schmidt等人的论文所述）
%-希尔伯特信封
%-基于小波变换的特征
%-基于PSD的特征
% This function was developed for use in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
%
%% INPUTS:
% audio_data: array of data from which to extract features
%从中提取特征的数据数组
% Fs: the sampling frequency of the audio data
%音频数据的采样频率
% figures (optional): boolean variable dictating the display of figures
%指示图形显示的布尔变量
%% OUTPUTS:
% PCG_Features: array of derived features
%派生特征数组
% featuresFs: the sampling frequency of the derived features. This is set
% in default_Springer_HSMM_options.m
%派生特征的采样频率。这是在默认的Springer_HSMM_选项中设置的
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

function [PCG_Features, featuresFs] = getSpringerPCGFeatures(audio_data, Fs, figures)
% function PCG_Features = getSpringerPCGFeatures(audio, Fs)
% Get the features used in the Springer segmentation algorithm
%获取Springer分割算法中使用的特征



if(nargin < 3)
    figures = 0;
end

springer_options = default_Springer_HSMM_options;


% Check to see if the Wavelet toolbox is available on the machine:
%检查机器上是否有小波工具箱：
include_wavelet = springer_options.include_wavelet_feature;
featuresFs = springer_options.audio_segmentation_Fs; % Downsampled feature sampling frequency
% 低通是50

%% 25-400Hz 4th order Butterworth band pass  25-400Hz四阶巴特沃斯带通
audio_data = butterworth_low_pass_filter(audio_data,2,400,Fs, false);
audio_data = butterworth_high_pass_filter(audio_data,2,25,Fs);

%% Spike removal from the original paper:   去除原文件的噪音
audio_data = schmidt_spike_removal(audio_data,Fs);



%% Find the homomorphic envelope            找到同态包络
% resample为信号降采样处理，理解如下：
% B=resample(x,90,250);  % 采样从250Hz降到90Hz，如果250在前,就是插值从90到250,可
% 以看B的长度,250Hz采样4000个数据等于90hz采样1440个数据,这就是降采样。这里是：
%心音时长X50/1000

homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs);
% Downsample the envelope:                       降低采样包络
downsampled_homomorphic_envelope = resample(homomorphic_envelope,featuresFs, Fs);
% normalise the envelope:                        正则化包络
downsampled_homomorphic_envelope = normalise_signal(downsampled_homomorphic_envelope);


%% Hilbert Envelope                         希尔伯特包络
hilbert_envelope = Hilbert_Envelope(audio_data, Fs);
downsampled_hilbert_envelope = resample(hilbert_envelope, featuresFs, Fs);
downsampled_hilbert_envelope = normalise_signal(downsampled_hilbert_envelope);

%% Power spectral density feature:         功率谱密度特征：

psd1 = get_PSD_feature_Springer_HMM(audio_data, Fs, 40,60)';
psd = resample(psd1, length(downsampled_homomorphic_envelope), length(psd1));
psd = normalise_signal(psd);

%% Wavelet features:                       小波特征：

if(include_wavelet)
    wavelet_level = 3;
    wavelet_name ='rbio3.9';
    
    % Audio needs to be longer than 1 second for getDWT to work:
    %音频需要超过1秒才能使getDWT工作
    if(length(audio_data)< Fs*1.025)
        audio_data = [audio_data; zeros(round(0.025*Fs),1)];
    end
    
    [cD, cA] = getDWT(audio_data,wavelet_level,wavelet_name);
    
    wavelet_feature = abs(cD(wavelet_level,:));
    wavelet_feature = wavelet_feature(1:length(homomorphic_envelope));
    downsampled_wavelet = resample(wavelet_feature, featuresFs, Fs);
    downsampled_wavelet =  normalise_signal(downsampled_wavelet)';

end


%% 经验模态分解
y = emd(audio_data,'fix',10);
EMD = resample(y, length(downsampled_homomorphic_envelope), length(y));
EMD_normalise= normalise_signal(EMD);
EMD_normalise = EMD_normalise';
EMD_normalise = resample(EMD_normalise, featuresFs, Fs);
EMD_normalise = normalise_signal(EMD_normalise);

%%

% if(include_wavelet)
    PCG_Features = [downsampled_homomorphic_envelope,downsampled_hilbert_envelope,downsampled_wavelet,psd,EMD_normalise];
% else
%     PCG_Features = [downsampled_homomorphic_envelope, downsampled_hilbert_envelope, psd];
% end

%% Plotting figures
% if(figures)
%     figure('Name', 'EMD');
% %     plot(audio_data);
% %     hold on;
%     plot(yyy,'r');
% %     legend('Original Signal',' wavelet feature')
if(figures)
    figure('Name', 'PCG features');
    t1 = (1:length(audio_data))./Fs;
    plot(t1,audio_data);
    hold on;
    t2 = (1:length(PCG_Features))./featuresFs;
    plot(t2,PCG_Features);
    legend('audio data','Homomorphic Envelope','Hilbert Envelope','PSD Feature','wavelet feature')
%     pause();
end