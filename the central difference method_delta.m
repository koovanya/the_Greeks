S0=48;
r=0.05;
sigma=0.3;
T=1;
n=100000;
K=52;
h=[1 0.1 0.01];

X=zeros(3,n);
Y=zeros(3,n);
H=zeros(3,n);
Z=zeros(3,n);
delta_FDM=zeros(1,3);
se=zeros(1,3);

d1=(log(S0./K)+(r+0.5*sigma.^2)*T)./(sigma*sqrt(T));
delta_BSM=normcdf(d1)-1;


for k=1:3
    for i=1:n
        Z(k,i)=randn;
        X(k,i)=(S0+h(k))*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*Z(k,i));
        Y(k,i)=(S0-h(k))*exp((r-0.5*sigma^2)*T+sigma*sqrt(T)*Z(k,i));
        H(k,i)=(exp(-r*T)*max(K-X(k,i),0)-exp(-r*T)*max(K-Y(k,i),0))/(2*h(k));
    end
    delta_FDM(k)=mean(sum(H));
    se(k)=sqrt(sum(H(k,:))^2-n*delta_FDM(k)^2)/(n*(n-1));
end

