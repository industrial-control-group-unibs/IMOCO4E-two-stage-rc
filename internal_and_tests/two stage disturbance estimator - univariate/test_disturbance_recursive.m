clear all;close all;clc;

tc=1e-3; % sample period
s=tf('s');
z=tf('z',tc);
% process definition
while true
    P=rss(30); % random 30-order process
    P=P*exp(-0.1*s); % add delay

    if not(any(real(zero(P))>0)) % the method does not work if there are zero with positive real part (inversion required)
        break
    end
end

Pd=(c2d(P,tc)); 

t=(0:tc:10)';


%%
std_n=0.01; % measure noise

n=std_n*randn(length(t),1);

u=square(2*pi*t/2);

v_=sin(2*pi*t/(0.2*t(end))); % variation of the independent variable 

%%
Ps=delay2z(Pd); % considere delay as z^-d

rel_order=length(pole(Ps))-length(zero(Ps));

F=1/(z^rel_order); % transfer function with delay
Pinv=Ps^-1*F; % inverse of the transer function without delay

%%

sigma_noise=std_n; % estimation of the noide

sigma=1;       % rbf parameters
l=1;           % rbf parameters
tau=[sigma;l]; % rbf parameters

nt=5;  % size of time-based GP
nv=50; % size of spatial GP

v_tbuf=zeros(nt,1);
d_tbuf=zeros(nt,1);

v_vbuf=linspace(-1,1,nv)';
d_vbuf=zeros(nv,1);
