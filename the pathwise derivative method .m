S0=48;
r=0.05;
sigma=0.3;
T=1;
n=100000;
K=48;
m=12;

S=zeros(m,n);
S(1,:)=S0;
S_bar=zeros(1,n);
delta=zeros(1,n);
dt=T/m;

for j=1:n
    for i=2:m
        S(i,j)=S(i-1,j)*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*randn);
    end
    S_bar(j)=sum(S(:,j))/m;
    delta(j)=exp(-r*T)*(S_bar(j)>K)*S(m,j)/S0;
end

delta_pathwise=mean(delta);
se=sqrt(var(delta)/n);
