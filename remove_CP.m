function y=remove_CP(x,Ncp,Noff)
% Remove cyclic prefix

%MIMO-OFDM Wireless Communications with MATLAB��   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

if nargin<3,  Noff=0;  end
y=x(:,Ncp+1-Noff:end-Noff);