%%%%%%%%%%%%%%%%%%%%%        ����CP��STO�����Ŷ�ʱƫ�    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        STO_estimation_sim1.m            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      date:2020��11��20��  author:����󽫾�   %%%%%%%%%%

%%%%%%%%%%%%%%%%%����˵��
%%%%����CP�����������غ���С�����㷨���STO�Ĺ���


%%%%%%    ���滷��
%����汾��MATLAB R2019a

clear, figure(1), clf, figure(2), clf
nSTOs = [-3 -3 2 2]; %��ӦSTO�Ĳ�����
CFOs = [0 0.5 0 0.5]; %CFO����
SNRdB = 30;   %SNR
MaxIter = 10;  %��������
%CFOs = [0 0 0 0];
Nfft = 128; %FFT��С
Ng = Nfft/4; % GI����
Nofdm = Nfft+Ng; % OFDM���ų���
Nbps = 2; % 2/4 ��Ӧ QPSK/16QAM 
M = 2^Nbps; 
Es = 1; 
A = sqrt(3/2/(M-1)*Es); % QAM��һ������
N = Nfft; 
com_delay = Nofdm/2;
Nsym = 100;
rand('seed',1); 
randn('seed',1);
for i=1:length(nSTOs)
   nSTO= nSTOs(i);  
   CFO= CFOs(i);
   x = []; % ��ʼ���źſ�
   for m=1:Nsym % random bits generates 
      msgint=randi([0 M-1],1,N); %bits_generator(1,Nsym*N,Nbps) 
      Xf = A.*qammod(msgint,M,'UnitAveragePower',true);
      xt = ifft(Xf,Nfft);  
      x_sym = add_CP(xt,Ng); %��CP
      x = [x x_sym]; 
   end
   %*********************** �ŵ� ************************%
   %%%%%�����������������ŵ����ȼ�����û���ŵ�
   y = x;  % û���ŵ�Ӱ��
   
   sig_pow = y*y'/length(y); % sig_pow= mean(mean(y.*conj(y),2))
   
   %%%%%Ƶ��ƫ�� +  ���Ŷ�ʱƫ�� 
   y_CFO= add_CFO(y,CFO,Nfft); %��CFO
   y_CFO_STO= add_STO(y_CFO,-nSTO);  %��STO
   v_ML=zeros(1,Ng); 
   v_Cl=zeros(1,Ng);
   Mag_cor= 0; 
   Mag_dif= 0;
   %%Add additive white gaussian noise
   for iter=1:MaxIter
      %%%%%%%%%%%������
      y_aw = awgn(y_CFO_STO,SNRdB,'measured');
      
      %%%%%%%Symbol Timing Acqusition
      [STO_cor,mag_cor]= STO_by_correlation(y_aw,Nfft,Ng,com_delay); %�����Դ�
      [STO_cor_temp,mag_cor_temp]= STO_by_correlation_sim1(y_aw,Nfft,Ng,com_delay); %���Լ���д
      %%%%%����֤���������ߺ������һ��
      
      
      [STO_dif,mag_dif] = STO_by_difference(y_aw,Nfft,Ng,com_delay); %�����Դ�
      [STO_dif_temp,mag_dif_temp] = STO_by_difference_sim1(y_aw,Nfft,Ng,com_delay); %���Լ���д
      %%%%%����֤���������ߺ������һ��
      
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
   legend('������ص�','���ڲ�ֵ��С��'); 
end % End of for loop of i