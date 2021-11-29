function CFO_est=CFO_Moose_sim1(y,Nfft,Nps)

CFO_est_temp = zeros(1,2*Nps+1);
for i=0:2*Nps
    D = Nfft/Nps;
    y_temp1 = y(D*i+1:D*i+D);
    y_temp2 = y(D*(i+1)+1:D*(i+1)+D);
    CFO_est_temp(i+1) = Nps*angle(y_temp2*y_temp1')/(2*pi);  
end
CFO_est = mean(CFO_est_temp);

