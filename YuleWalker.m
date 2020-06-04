%%  Yule-Walker for ARMA model
% Author: Xianrui Wang, 2019100467
% Center of Intelligent Acoustics and Immersive Communications 
% Contact me: wangxianrui@mail.nwpu.edu.cn
%-------------------------------------------------------------------------
% smapleNums: number of sampling points
% modelOrder: order of AR model
% R: auto correlation of signal
% Rx: the second half of R and have legnth of modelOrder
% A: Yule Walker equation can be denoted by p=Aa, a is the parameter need 
% to be estimated, and [A]ij=R(i-j),[p]i=R(i)
% H: system transfer function = G/(1+\sum_{i=1}^{modelOrder}a(i)z^{-i})
% G: system gain = Rx(0)+\sum_{i=1}^{modelOrder}a(i)Rx(i)
% a: coefficients of AR model using direct inverse
% b: coefficients of AR model using Levenson algorithm
%------------------------------------------------------------------------
clear;clc;
%% Basic parameter initial 
smapleNums = 128;
modelOrder = 64;
t = 0:smapleNums-1;
inputSnr = 30;
x =10*sin(2*pi*0.1*t+pi/3) + 10*sin(2*pi*0.2*t+pi/4);
xn=awgn(x,inputSnr);
R = xcorr(xn,'biased');
Rx = zeros(1,modelOrder+1);
for i = 1:modelOrder+1
    Rx(i) = R(i+smapleNums-1);
end
%-------------------------------------------------------------------------
%% Yule Walker equation
A = zeros(modelOrder,modelOrder);
for m = 1:modelOrder
    for n = 1:modelOrder
        A(m,n) = Rx(max(m,n)-min(m,n)+1);
    end
end
p = zeros(modelOrder,1);
for m = 1:modelOrder
    p(m,1) = -Rx(m+1);
end
%% Direct inverse solution of Yule walker equation 
a = (A\p).';
G_squre = Rx(1);
for i = 1:modelOrder      
    G_squre = G_squre+a(i)*Rx(i);
end
G=G_squre^0.5;
fprintf('System Gain G=%f\n',G);
[H,w] = freqz(G,[1,a],smapleNums);
Hf = abs(H);
%--------------------------------------------------------------------------
%% Solution with Levinson-Durbin algorithm
% this part comes directly from slides of chapter 5
r=Rx.';
b=1;
epsilon=r(1);
for j=2:modelOrder+1
    gamma=-r(2:j).'*flipud(b)/epsilon;
    b=[b;0]+gamma*[0;conj(flipud(b))];
    epsilon = epsilon*(1-abs(gamma)^2);
end
b=b(2:65);
[H_Levinson,w] = freqz(1,[1;b],smapleNums);
Hf_Levinson = abs(H_Levinson);
%% Plot transfer function
f_range = w/(2*pi);
figure 
plot(f_range,Hf,'LineWidth',2);
xlabel('Frequency (kHz)');
title('Direct Inverse solution');
grid on;
figure 
plot(f_range,Hf_Levinson,'LineWidth',2);
xlabel('Frequency (kHz)');
title('Levinson-Durbin algorithm');
grid on;
