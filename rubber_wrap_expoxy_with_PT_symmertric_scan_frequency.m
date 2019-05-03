clear all;
clc;
% parametric input
K=1.25;
M=2.0;
m=0.5;
alpha=0.5;
G=0.3;
d=0.1;
mag=[];
magf=[];
kw=[];
fr=[];
%resonant frequency
w_0_loss = sqrt(2*G*(1+alpha*i)/m);
w_0_gain = sqrt(2*G*(1-alpha*i)/m);
%扫描角频率
w=0:0.001:2;
%虚部不为0的频率集合
jg=0;
%虚部为0的频率集合
jf=0;
for h=1:length(w)
    %dispersion relation
    q = acos((2*K-w(h)^2*(M+m*w_0_loss^2/(w_0_loss^2-w(h)^2)))*(2*K-w(h)^2*(M+m*w_0_gain^2/(w_0_gain^2-w(h)^2)))/(2*K^2)-1);
%     if imag(q) ~= 0 
%         jg=jg+1;
%         mag(jg)=imag(q);
%         magf(jg)=w(h);
%     else
%         jf=jf+1;
%         kw(jf)=q;
%         fr(jf)=w(h);
%     end

%全部波矢的实部和虚部分别取出并画图
jg=jg+1;
mag(jg)=imag(q);
magf(jg)=w(h);
jf=jf+1;
kw(jf)=real(q)/(pi);
fr(jf)=w(h);
end
%画图，色散曲线
figure;
plot(kw,fr,'LineWidth',2);
hold on;
plot(mag,magf,'r','LineWidth',2);