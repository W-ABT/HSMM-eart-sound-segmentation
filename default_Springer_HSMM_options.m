% function springer_options = default_Springer_HSMM_options()
%
% The default options to be used with the Springer segmentation algorithm.
% USAGE: springer_options = default_Springer_HSMM_options
%
%��˹����ָ��㷨һ��ʹ�õ�Ĭ��ѡ��÷���springer_options = default_Springer_HSMM_options
%
% Developed for use in the paper:
% D. Springer et al., "Logistic Regression-HSMM-based Heart Sound 
% Segmentation," IEEE Trans. Biomed. Eng., In Press, 2015.
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

function springer_options = default_Springer_HSMM_options()

%% The sampling frequency at which to extract signal features:
%��ȡ�ź������Ĳ���Ƶ�ʣ�
springer_options.audio_Fs = 1000;

%% The downsampled frequency
%Set to 50 in Springer paper
%�²���Ƶ�ʣ���springerֽ��������Ϊ50
springer_options.audio_segmentation_Fs = 50;


%% Tolerance for S1 and S2 localization
%s1��s2��λ����
springer_options.segmentation_tolerance = 0.1;%seconds

%% Whether to use the mex code or not:
% The mex code currently has a bug. This will be fixed asap.
%�Ƿ�ʹ��MEX���룺
%MEX���뵱ǰ��һ�������������ᾡ�����ġ�
springer_options.use_mex = true;

%% Whether to use the wavelet function or not:
%�Ƿ�ʹ��С��������
springer_options.include_wavelet_feature = true;

