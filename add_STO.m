function y_STO=add_STO(y,iSTO)

%MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

if iSTO>=0
    y_STO = [y(iSTO+1:end) zeros(1,iSTO)]; % advance
else
    y_STO = [zeros(1,-iSTO) y(1:end+iSTO)];  % delay
end