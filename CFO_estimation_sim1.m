%%%%%%%%%%%%%%%%%%%%%        ����CP��CFO���ز�Ƶƫ���ƣ�    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%        CFO_estimation_sim1.m            %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      date:2020��12��14��  author:����󽫾�   %%%%%%%%%%

%%%%%%%%%%%%%%%%%����˵��
%%%���ʱ�����CP�ķ�����Ƶ���Moose/Classen���������ں���CFO����


%%%%%%    ���滷��
%����汾��MATLAB R2019a
clear, clf
CFO = 0.15;
% CFO = 0;
Nfft=128; % FFT size   
Nbps=2; 
M=2^Nbps; % Number of bits per (modulated) symbol
Es=1; 
A=sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor
N=Nfft;
Ng=Nfft/4;
Nofdm=Nfft+Ng; 
Nsym=3;  
% h=complex(randn,randn)/sqrt(2);
% %h=[1 zeros(1,5)]; 
% channel(h,0);  
%Transmit signal
x=[];
for m=1:Nsym
   msgint=randi([0 M-1],1,N); %bits_generator(1,Nsym*N,Nbps) 

   if m<=2
       Xp = add_pilot(zeros(1,Nfft),Nfft,4); 
       Xf=Xp; % add_pilot
   else  %Xf= QAM(msgint((i-1)*N+1:i*N),Nbps);  % constellation mapping. average power=1        
       Xf = A.*qammod(msgint,M,'UnitAveragePower',true);
   end                        
   xt = ifft(Xf,Nfft);  
   x_sym = add_CP(xt,Ng);
   x= [x x_sym];
end    

%*************************  �ŵ� **************
%channel ����������ŵ�
y=x; % No channel effect

%Signal power calculation
sig_pow= y*y'/length(y); % Signal power calculation

%%%%
SNRdBs= 0:3:30;  
% SNRdBs= 100; ��100��Ϊ���Գ���  
MaxIter = 100;  
for i=1:length(SNRdBs)
   SNRdB = SNRdBs(i);
   MSE_CFO_CP = 0; 
   MSE_CFO_Moose = 0; 
   MSE_CFO_Classen = 0;
   rand('seed',1); 
   randn('seed',1);  % Initialize seed for random number generator
   y_CFO= add_CFO(y,CFO,Nfft); % Add CFO
   %%%%��ε���ȡƽ��
   for iter=1:MaxIter
      %y_aw=add_AWGN(y_CFO,sig_pow,SNRdB,'SNR',Nbps);  % AWGN added, signal power=1
      y_aw = awgn(y_CFO,SNRdB,'measured');  % AWGN added, signal power=1
      
      Est_CFO_CP = CFO_CP(y_aw,Nfft,Ng); % CP-based 
      MSE_CFO_CP = MSE_CFO_CP + (Est_CFO_CP-CFO)^2;
      
      Est_CFO_Moose = CFO_Moose(y_aw,Nfft); % Moose (based on two consecutive preambles)
      MSE_CFO_Moose = MSE_CFO_Moose + (Est_CFO_Moose-CFO)^2;
      
      Est_CFO_Classen = CFO_Classen(y_aw,Nfft,Ng,Xp); % Classen (Pilot-based)
      MSE_CFO_Classen = MSE_CFO_Classen + (Est_CFO_Classen-CFO)^2;
      
   end % the end of for (iter) loop
   MSE_CP(i) = MSE_CFO_CP/MaxIter; 
   MSE_Moose(i) = MSE_CFO_Moose/MaxIter;  
   MSE_Classen(i) = MSE_CFO_Classen/MaxIter;
end%ebn0 end    
semilogy(SNRdBs, MSE_CP,'-+');
grid on, hold on
semilogy(SNRdBs, MSE_Moose,'-x'); semilogy(SNRdBs, MSE_Classen,'-*');
xlabel('SNR[dB]'), ylabel('MSE'); title('CFO Estimation'); %axis([0 30 10e-8 10e-2])
% str=sprintf('CFO = %1.2f',CFO);
legend('CP-based technique','Moose (Preamble-based)','Classen (Pilot-based)');
% legend(str);

%%********************* ʵ���� *************
%%%��¼��CFO_fig1���Ա��˲�ͬ�㷨��CFO���ƽ��
%%%2020��12��14�� ����������ȷ
%%%�ɹ�����CFOС�����֣�CFO���ƺ����ٿ���
