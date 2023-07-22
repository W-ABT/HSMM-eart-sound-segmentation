% function high_pass_filtered_signal = butterworth_high_pass_filter(original_signal,order,cutoff,sampling_frequency)
%
% High-pass filter a given signal using a forward-backward, zero-phase
% butterworth filter.
%高通滤波器使用正反向、零相位巴特沃斯滤波器对给定信号进行滤波。
%
%% INPUTS:
% original_signal: The 1D signal to be filtered                            原始信号：要滤波的1D信号
% order: The order of the filter (1,2,3,4 etc). NOTE: This order is        顺序：过滤器的顺序（1、2、3、4等）。
% effectively doubled as this function uses a forward-backward filter that 注：由于此功能使用了确保零相位失真的
% ensures zero phase distortion                                            前向-后向滤波器，因此该阶数有效地加倍
% cutoff: The frequency cutoff for the high-pass filter (in Hz)            截止：高通滤波器的频率截止值（Hz）
% sampling_frequency: The sampling frequency of the signal being filtered  采样频率：被滤波信号的采样频率（Hz）。
% (in Hz).
% figures (optional): boolean variable dictating the display of figures    图形（可选）：指示图形显示的布尔变量
%
%% OUTPUTS:
% high_pass_filtered_signal: the high-pass filtered signal.                高通滤波信号：高通滤波信号。
%
% This code is derived from the paper:
% S. E. Schmidt et al., "Segmentation of heart sound recordings by a
% duration-dependent hidden Markov model," Physiol. Meas., vol. 31,
% no. 4, pp. 513-29, Apr. 2010.
%
% Developed by David Springer for comparison purposes in the paper:
% D. Springer et al., ?Logistic Regression-HSMM-based Heart Sound
% Segmentation,? IEEE Trans. Biomed. Eng., In Press, 2015.
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

function high_pass_filtered_signal = butterworth_high_pass_filter(original_signal,order,cutoff,sampling_frequency, figures)

if nargin < 5
    figures = 0;
end

%Get the butterworth filter coefficients
[B_high,A_high] = butter(order,2*cutoff/sampling_frequency,'high');

%Forward-backward filter the original signal using the butterworth
%coefficients, ensuring zero phase distortion
high_pass_filtered_signal = filtfilt(B_high,A_high,original_signal);

if(figures)
    
    figure('Name','High-pass filter frequency response');
    [sos,g] = zp2sos(B_high,A_high,1);	     % Convert to SOS form
    Hd = dfilt.df2tsos(sos,g);   % Create a dfilt object
    h = fvtool(Hd);	             % Plot magnitude response
    set(h,'Analysis','freq')	     % Display frequency response
    
    figure('Name','Original vs. high-pass filtered signal');
    plot(original_signal);
    hold on;
    plot(high_pass_filtered_signal,'r');
    legend('Original Signal', 'High-pass filtered signal');
    pause();
end

