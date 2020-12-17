function CFO_est=CFO_Moose(y,Nfft)
% Frequency-domain CFO estimation using Moose method based on two consecutive identical OFDM symbols 

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

for i=0:1   
    
    %%自己改动
%  
%     y_temp = y( (Nfft+32)*i + 33 : (Nfft+32)*i+ 160); 
%    Y(i+1,:)= fft(y_temp,Nfft);

% %%书中自带公式
    y_temp = y(Nfft*i+1:Nfft*(i+1));
    Y(i+1,:)= fft(y_temp,Nfft);
    
% %%%
%     Y(i+1,:)= fft(y(Nfft*i+1 : Nfft*(i+1)+31),Nfft);   
end
CFO_est = angle(Y(2,:)*Y(1,:)')/(2*pi); % Eq.(5.30)
%fprintf(' FFO-Esti_FFO=%8.5f-%8.5f\n', FFO,Esti_FFO);
% Esti_FFO=Esti_FFO*(Nfft/(Nfft+Nfft/4));
%FFO_MSE = FFO_MSE + (Esti_FFO-FFO)*(Esti_FFO-FFO);
