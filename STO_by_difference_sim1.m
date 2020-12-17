function [STO_est,Mag] = STO_by_difference_sim1(y,Nfft,Ng,com_delay)
% STO estimation by minimizing the difference between CP and rear part of OFDM symbol
% estimates STO by minimizing the difference between CP (cyclic prefix) 
%     and rear part of OFDM symbol
% Input:  y          = Received OFDM signal including CP
%          Ng         = Number of samples in CP (Guard Interval)
%          com_delay = Common delay
% Output: STO_est   = STO estimate
%           Mag        = Correlation function trajectory varying with time

N_ofdm = Nfft+Ng; 
minimum = 100; 
STO_est = 0;
if nargin<4
    com_delay = N_ofdm/2; 
end

temp = abs(y(com_delay : com_delay + Ng-1)) - abs(y(com_delay+Nfft : com_delay + Ng-1 +Nfft));
minimum = temp*temp'; % ¹«Ê½.(5.12)

for k = 1:N_ofdm
    temp = abs(y(com_delay + k  : com_delay + k + Ng-1)) - abs(y(com_delay+Nfft+k : com_delay + Ng-1 +Nfft+k));
    Mag(k) = temp*temp';
    if Mag(k)<minimum
        minimum =  Mag(k);
        STO_est = N_ofdm - com_delay - k + 1;
    end
end
end

