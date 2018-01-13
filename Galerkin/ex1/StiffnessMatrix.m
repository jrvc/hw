function K=StiffnessMatrix(n)
%calculation of the stiffness matrix K
K=zeros(n+1,n+1);
K(1,1)=2;
for i=1:1:n
    for j=1:1:n 
        K(i+1,j+1)=2^(i+j-1)*(i*j/(i+j-1)+4/(i+j+1));
    end
    K(i+1,1)=2^(i+1)/(i+1);
    K(1,i+1)=2^(i+1)/(i+1);
end
