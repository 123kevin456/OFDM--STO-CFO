%%%%%%%%%%%%%%%%%%%%%        基于CP的STO（符号定时偏差）    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        STO_estimation_sim1.m            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      date:2020年11月20日  author:飞蓬大将军   %%%%%%%%%%

%%%%%%%%%%%%%%%%%程序说明
%%%%基于CP，采用最大相关和最小距离算法完成STO的估计


%%%%%%    仿真环境
%软件版本：MATLAB R2019a

clear, figure(1), clf, figure(2), clf
nSTOs = [-3 -3 2 2]; %对应STO的采样数
CFOs = [0 0.5 0 0.5]; %CFO向量
SNRdB = 30;   %SNR
MaxIter = 10;  %迭代次数
%CFOs = [0 0 0 0];
Nfft = 128; %FFT大小
Ng = Nfft/4; % GI长度
Nofdm = Nfft+Ng; % OFDM符号长度
Nbps = 2; % 2/4 对应 QPSK/16QAM 
M = 2^Nbps; 
Es = 1; 
A = sqrt(3/2/(M-1)*Es); % QAM归一化因子
N = Nfft; 
com_delay = Nofdm/2;
Nsym = 100;
rand('seed',1); 
randn('seed',1);
for i=1:length(nSTOs)
   nSTO= nSTOs(i);  
   CFO= CFOs(i);
   x = []; % 初始化信号块
   for m=1:Nsym % random bits generates 
      msgint=randi([0 M-1],1,N); %bits_generator(1,Nsym*N,Nbps) 
      Xf = A.*qammod(msgint,M,'UnitAveragePower',true);
      xt = ifft(Xf,Nfft);  
      x_sym = add_CP(xt,Ng); %加CP
      x = [x x_sym]; 
   end
   %*********************** 信道 ************************%
   %%%%%在这里根据需求添加信道，先假设是没有信道
   y = x;  % 没有信道影响
   
   sig_pow = y*y'/length(y); % sig_pow= mean(mean(y.*conj(y),2))
   
   %%%%%频率偏移 +  符号定时偏移 
   y_CFO= add_CFO(y,CFO,Nfft); %加CFO
   y_CFO_STO= add_STO(y_CFO,-nSTO);  %加STO
   v_ML=zeros(1,Ng); 
   v_Cl=zeros(1,Ng);
   Mag_cor= 0; 
   Mag_dif= 0;
   %%Add additive white gaussian noise
   for iter=1:MaxIter
      %%%%%%%%%%%加噪声
      y_aw = awgn(y_CFO_STO,SNRdB,'measured');
      
      %%%%%%%Symbol Timing Acqusition
      [STO_cor,mag_cor]= STO_by_correlation(y_aw,Nfft,Ng,com_delay); %书中自带
      [STO_cor_temp,mag_cor_temp]= STO_by_correlation_sim1(y_aw,Nfft,Ng,com_delay); %我自己编写
      %%%%%经验证，以上两者函数结果一致
      
      
      [STO_dif,mag_dif] = STO_by_difference(y_aw,Nfft,Ng,com_delay); %书中自带
      [STO_dif_temp,mag_dif_temp] = STO_by_difference_sim1(y_aw,Nfft,Ng,com_delay); %我自己编写
      %%%%%经验证，以上两者函数结果一致
      
      v_ML(-STO_cor+Ng/2)= v_ML(-STO_cor+Ng/2)+1;
      v_Cl(-STO_dif+Ng/2)= v_Cl(-STO_dif+Ng/2)+1;
      Mag_cor= Mag_cor + mag_cor; 
      Mag_dif= Mag_dif + mag_dif;
   end % End of for loop of iter
   %%%%%%% Probability
   v_ML_v_Cl= [v_ML; v_Cl]*(100/MaxIter);
   figure(1+i-1); 
   set(gca,'fontsize',9);
%    subplot(220+i)
   bar(-Ng/2+1:Ng/2,v_ML_v_Cl');
   hold on, grid on
   str = sprintf('nSTO Estimation: nSTO=%d, CFO=%1.2f, SNR=%3d[dB]',nSTO,CFO,SNRdB);           
   title(str); 
   xlabel('Sample'), ylabel('Probability');
   legend('ML','Classen'); 
   axis([-Ng/2-1 Ng/2+1 0 100])
   %%%%%%% Time metric
   Mag_cor = Mag_cor/MaxIter; 
   [Mag_cor_max,ind_max] = max(Mag_cor); 
   nc= ind_max-1-com_delay; 
   Mag_dif = Mag_dif/MaxIter; 
   [Mag_dif_min,ind_min] = min(Mag_dif); 
   nd= ind_min-1-com_delay;
   nn= -Nofdm/2 + [0:length(Mag_cor)-1]; 
   nt= nSTO;
%    figure(2);
%    subplot(220+i);
   figure(5+i-1); 
   plot(nn,Mag_cor,nn,1.5*Mag_dif,'r:','markersize',1);
   hold on
     
   stem(nc,Mag_cor_max,'b','markersize',5);
   stem(nSTO,Mag_cor(nSTO+com_delay+1),'k.','markersize',5); % Estimated/True Maximum value
   str1 = sprintf('STO Estimation - ML(b-)/Classen(r:) for nSTO=%d, CFO=%1.2f',nSTO,CFO); %,SNRdB);
   title(str1); 
   xlabel('Sample'), ylabel('Magnitude'); 
 
   %stem(n1,Mag_dif_min,'r','markersize',5)
   stem(nd,Mag_dif(nd+com_delay+1),'r','markersize',5);
   stem(nSTO,Mag_dif(nSTO+com_delay+1),'k.','markersize',5); % Estimated/True Minimum value
   set(gca,'fontsize',9, 'XLim',[-32 32], 'XTick',[-10 -3 0 2 10]); %, xlim([-50 50]),
   legend('基于相关的','基于差值最小的'); 
end % End of for loop of i