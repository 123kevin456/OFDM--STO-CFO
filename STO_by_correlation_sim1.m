%%%%%%%%%%%%%%%%%%%%%        STO相关函数    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        STO_by_correlation_sim1.m    %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      date:2020年11月20日  author:飞蓬大将军   %%%%%%%%%%


%%%%%%    仿真环境
%软件版本：MATLAB R2019a

function [STO_est, Mag]=STO_by_correlation_sim1(y,Nfft,Ng,com_delay)


% STO estimation by maximizing the correlation between CP and rear part of OFDM symbol
% estimates STO by maximizing the correlation between CP (cyclic prefix)  
%     and rear part of OFDM symbol
% Input:  y         = Received OFDM signal including CP
%         Ng        = Number of samples in Guard Interval (CP)
%         com_delay = Common delay
% Output: STO_est   = STO estimate
%         Mag       = Correlation function trajectory varying with time


N_ofdm = Nfft+Ng; 
if nargin<4
    com_delay = N_ofdm/2;
end

yy = y(com_delay : com_delay+ Ng-1)*y(Nfft+com_delay : Nfft+com_delay+Ng-1)';
maximum = abs(yy);

for k = 1:N_ofdm
    yy_temp = y(com_delay+k : com_delay+k+Ng-1)*y( Nfft+com_delay+k : Nfft+com_delay+k+Ng-1)'; %% 公式.(5.13)
    Mag(k) = abs(yy_temp);
    if abs(yy_temp) > maximum
        maximum = abs(yy_temp);
        STO_est = N_ofdm - com_delay - k + 1;
    end
end
end
