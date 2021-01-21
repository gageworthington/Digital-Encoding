ibs = [1,1,0,0,1,0,1,0]
tlc = 'nrz-s'
Rb = 1
A = 1

lc_gen_plot(ibs,tlc,Rb,A)

function [x T] = lc_gen_plot(ibs,tlc,Rb,A)
% This function implements various Line Coding techniques
% To use the code, change the parameters ibs, tlc, Rb, and A 
% to edit the input bistream, encoding type, bit rate, and 
% signal amplitude, respectively. 
%
% INPUT:
% tlc = string accepting the tlc of coding as the following:
% 'NRZ-L'
% 'nrz-m'
% 'nrz-s'
% 'RZ'
% 'manch'
% '__ami'
% 'bi_rz'
% 'd_man'
% ibs = input bits row vector
% Rb = Bit time rate
% A = Amplitude of the coding
% OUTPUT:
% x = Output line coding row vector
% T = time vector

clf
if nargin<4, A = 1;end
if nargin<3, Tb = 1e-9;end
if nargin<2, ibs = [1 0 1 0];end
if nargin<1, tlc = 'uninrz';end

%---Implementation starts here
Tb = 1/Rb; %---Bit period
Fs = 4*Rb;
N = length(ibs); %---Bit Length of input bits
tTb = linspace(0,Tb); %---interval of bit time period
x = [];
i = 0; %---Index for __ami
prevstate = 0; %---Index array for nrz-m
tlc = lower(tlc)
switch lower(tlc)

case 'nrz-l'
for k = 1:N
x = [x A*ibs(k)*ones(1,length(tTb))];
end

case 'nrz-m'
for k = 1:N
if k == 1
if ibs(k) == 1
x = [x A*ibs(k)*ones(1,length(tTb))];
else
x = [x A*ibs(k)*0*ones(1,length(tTb))];
end
continue
end
if ibs(k) == 1
x = [x A*ibs(k)*~x(100*(k-1))*ones(1,length(tTb))];
else
x = [x A*x(100*(k-1))*ones(1,length(tTb))];
end
end

case 'nrz-s'
for k = 1:N
if k == 1
if ibs(k) == 0
x = [x A*~ibs(k)*ones(1,length(tTb))];
else
x = [x A*ibs(k)*0*ones(1,length(tTb))];
end
continue
end
if ibs(k) == 0
x = [x A*~ibs(k)*~x(100*(k-1))*ones(1,length(tTb))];
else
x = [x A*x(100*(k-1))*ones(1,length(tTb))];
end
end

case '__ami'
for k = 1:N
if ibs(k) == 1
x = [x A*ibs(k)*((-1)^i)*ones(1,length(tTb))];
i = i + 1
else
x = [x A*ibs(k)*ones(1,length(tTb))];
end
end

case 'rz'
for k = 1:N
x = [x A*ibs(k)*ones(1,length(tTb)/2) 0*ibs(k)*ones(1,length(tTb)/2)];
end

case 'bi_rz'
for k = 1:N
if ibs(k) == 1
x = [x A*ibs(k)*((-1)^i)*ones(1,length(tTb)/2)];
x = [x A*ibs(k)*0*ones(1,length(tTb)/2)]
i = i + 1
else
x = [x A*ibs(k)*ones(1,length(tTb))];
end
end

case 'manch'
for k = 1:N
c = ones(1,length(tTb)/2);
b = -1*ones(1,length(tTb)/2);
p = [c b];
x = [x (-(-1)^(ibs(k)+1))*A/2*p];
end

case 'd_man'
for k = 1:N
if k == 1
% if k > 0 ---- use as a last resort
if ibs(k) == 0
x = [x A/2*ones(1,length(tTb)/2)]
x = [x -A/2*ones(1,length(tTb)/2)]
continue
else
x = [x -A/2*ones(1,length(tTb)/2)]
x = [x A/2*ones(1,length(tTb)/2)]
continue
end
end
if ibs(k) == 0
x = [x -1*x(100*(k-1))*ones(1,length(tTb)/2)]
x = [x x(100*(k-1))*ones(1,length(tTb)/2)]
else
x = [x x(100*(k-1))*ones(1,length(tTb)/2)]
x = [x -1*x(100*(k-1))*ones(1,length(tTb)/2)]
end
end
end

T = linspace(0,N*Tb,length(x)); %---Time vector for n bits
plot(T,x)
xlim([0-(0.05*Tb*N) Tb*N+(0.05*Tb*N)]);
if min(x) == 0
ylim([min(x)-0.05 max(x)+(0.05*max(x))])
else
ylim([min(x)-(abs(0.05*min(x))) max(x)+(0.05*max(x))])
end
% axis([0-(0.05*Tb*N) Tb*N+(0.05*Tb*N) -1.5 1.5]);
title(tlc)
xlabel("time (sec)")
ylabel("amplitude")
grid on
end 
