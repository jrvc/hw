function bh=RightHandSide(n)
bh=zeros(n+1,1);
for i=1:1:n
   % syms x
    %bh(i+1)=eval(-8*pi*i*int(sin(8*pi*x)*x^(i-1),x,[0,2])+int(cos(8*pi*x)*x^i,x,[0,2]));
    fun=@(x)cos(8*pi.*x).*x.^i-8*pi*i*sin(8*pi.*x).*x.^(i-1);
    bh(i+1)=integral(fun,0,2);
end
% syms x
% bh(1)=eval(int(cos(8*pi*x),x,[0,2]));
bh(1)=integral(@(x)cos(8*pi.*x),0,2);