% function [cD cA] = getDWT(X,N,Name)
%
% finds the discrete wavelet transform at level N for signal X using the
% wavelet specified by Name.
%使用名称指定的小波查找信号X的N级离散小波变换。
%
%% Inputs:
% X: the original signal
% N: the decomposition level
% Name: the wavelet name to use
%
% X： 原始信号
% N： 分解层次
% 名称：小波的名字

%% Outputs:
% cD is a N-row matrix containing the detail coefficients up to N levels
% cA is the same for the approximations
% 
% cD是一个N行矩阵，包含N个级别的细节系数
% 对于近似值，cA是相同的
% This code was developed by David Springer for comparison purposes in the
% paper:
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

function [cD cA] = getDWT(X,N,Name)


% No DWT available for Morlet - therefore perform CWT:
% Morlet无DWT可用-因此执行CWT：
if(strcmp(Name,'morl'))
    
    c = cwt(X,1:N,'morl');
    
    cD = c;
    cA = c;
else
    % Preform wavelet decomposition
    % 预小波分解
    [c,l] = wavedec(X,N,Name);
    
    % Reorder the details based on the structure of the wavelet decomposition (see help in wavedec.m)
    % 根据小波分解的结构对细节重新排序（参见wavedec.m中的帮助）
    len = length(X);
    cD = zeros(N,len);
    for k = 1:N
        d = detcoef(c,l,k);
        d = d(:)';
        d = d(ones(1,2^k),:);
        cD(k,:) = wkeep1(d(:)',len);
    end
    cD =  cD(:);
    
    % Space cD according to spacing of floating point numbers:
    %  根据浮点数的间距将cD隔开： 
    I = find(abs(cD)<sqrt(eps));
    cD(I) = zeros(size(I));
    cD = reshape(cD,N,len);
    % cD = wcodemat(cfd,nbcol,'row');
    
    
    % Reorder the approximations based on the structure of the wavelet decomposition (see help in wavedec.m)
    % 根据小波分解的结构对近似值重新排序（参见wavedec.m中的帮助）
    len = length(X);
    cA = zeros(N,len);
    for k = 1:N
        a = appcoef(c,l,Name,k);
        a = a(:)';
        a = a(ones(1,2^k),:);
        cA(k,:) = wkeep1(a(:)',len);
    end
    cA =  cA(:);
    I = find(abs(cA)<sqrt(eps));
    cA(I) = zeros(size(I));
    cA = reshape(cA,N,len);
end