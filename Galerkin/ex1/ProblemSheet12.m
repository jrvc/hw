clear all;
clc;
k=8;
N=2^k;
h=2/N;
X=0:h:2;
tmp=ones(N+1,1);
set(0,'RecursionLimit',10*N)

phi0=@(x)((x>=X(1)).*(x<=X(2))*1/h.*(X(2)-x));
phiN=@(x)((x>=X(N)).*(x<=X(N+1))*1/h.*(x-X(N)));
phi=@(x,i)((x>=X(i-1)).*(x<X(i))*1/h.*(x-X(i-1))+(x>=X(i)).*(x<=X(i+1))*1/h.*(X(i+1)-x));

K=diag((2/h+2/3*h).*tmp,0)+diag((-1/h+h/6).*tmp(1:N,1),1)+diag((-1/h+h/6).*tmp(1:N,1),-1);
K(1,1)=K(1,1)/2;
K(N+1,N+1)=K(N+1,N+1)/2;



% K(1,2)=integral(@(x)(phi0(x).*phi(x,2)),X(1),X(2));
% K(2,1)=K(1,2);
% K(N+1,N)=integral(@(x)(phiN(x).*phi(x,N)),X(N),X(N+1));
% K(N,N+1)=K(N+1,N);

cond(K);


b=zeros(N+1,1);
b(1)=integral(@(x)(1+64*pi*pi).*cos(8*pi.*x).*phi0(x),X(1),X(2));
b(N+1)=integral(@(x)(1+64*pi*pi).*cos(8*pi.*x).*phiN(x),X(N),X(N+1));
for i=2:N
   b(i)=integral(@(x)(1+64*pi*pi).*cos(8*pi.*x).*phi(x,i),X(i-1),X(i+1));
end

Uh=K\b;

U=@(x)(Uh(1).*phi0(x)+Uh(N+1).*phiN(x));
for i=2:N
    U=@(x)(U(x)+Uh(i).*phi(x,i));
end
[A,B]=fplot(U,[0,2]);
V=cos(8*pi.*A);
plot(A,V,A,B)


