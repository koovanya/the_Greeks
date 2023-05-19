S0=48;
r=0.05;
sigma=0.3;
T=1;
n=100000;
K=48;

ST=zeros(1,n);
delta=zeros(1,n);
LR=zeros(1,n);

for i=1:n
    ST(i)=S0*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*randn);
    LR(i)=(log(ST(i)/S0)-(r-0.5*sigma^2)*T)/(S0*sigma^2*T);
    delta(i)=(ST(i)>K)*LR(i);
end

delta_LRM=mean(delta);
se=sqrt(var(delta)/n);
