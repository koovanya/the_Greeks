S0=48;
r=0.05;
sigma=0.3;
T=1;
n=100000;
K=52;
h=[0.1 0.01 0.001];

X=zeros(3,n);
Y=zeros(3,n);
H=zeros(3,n);
Z=zeros(3,n);
vega_FDM=zeros(1,3);
se=zeros(1,3);

d1=(log(S0./K)+(r+0.5*sigma.^2)*T)./(sigma*sqrt(T));
vega_BSM=S0*sqrt(T)*normpdf(d1);


for k=1:3
    for i=1:n
        Z(k,i)=randn;
        X(k,i)=S0*exp((r-0.5*(sigma+h(k))^2)*T+(sigma+h(k))*sqrt(T)*Z(k,i));
        Y(k,i)=S0*exp((r-0.5*(sigma-h(k))^2)*T+(sigma-h(k))*sqrt(T)*Z(k,i));
        H(k,i)=(exp(-r*T)*max(X(k,i)-K,0)-exp(-r*T)*max(Y(k,i)-K,0))/(2*h(k));
    end
    vega_FDM(k)=mean(sum(H));
    se(k)=sqrt(sum(H(k,:))^2-n*vega_FDM(k)^2)/(n*(n-1));
end

