%the function computes the approximate solution Uh of the problem
% a(u,phi)=(f,phi) for all phi belonging to H1([0,2])
%the discretezation of the space is given by Vh=span(1,x^2, ... x^n)
k=33;
% l2=zeros(1,k);
% for i=k:1:k
%    Ki=StiffnessMatrix(i);
%    K=Ki;
%    l2(i)=cond(Ki);
% end
% x=[1:1:k];
% plot(x,l2)
% title('condition number l2 for the stiffness matrices')
% xlabel('n')
% ylabel('l2(n)')
% grid on
K=StiffnessMatrix(k);
X=[0:0.001:2];
U=cos(8*pi.*X);

for i=k:1:k
    Ui=linsolve(K(1:i+1,1:i+1),RightHandSide(i));
    Yi=zeros(size(X));
    for j=1:1:i+1
        Yi=Yi+Ui(j).*X.^(j-1);
    end
    figure(i+1)
    plot(X,U,X,Yi)
end



