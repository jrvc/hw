close all;
clear all;
clc;

L = 1;
Mx = 100;
h = L / Mx;

T = 50;
Nt = 10000;
k = T / Nt;

lambda = k / (h*h);
theta = 1 / 2;

x = [0:h:L]';
g = sin(pi*x.*x);

A = zeros(Mx-1, Mx-1);

nOnes = ones(Mx-1, 1);
A = diag(-1-2*theta*lambda * nOnes, 0) + diag(theta*lambda*nOnes(1:Mx-1-1), -1) + diag(theta*lambda*nOnes(1:Mx-1-1), 1);
b = zeros(Mx-1,1);
u = g(2:Mx);
A = sparse(A);

for i = 1:Mx-1
    b(i) =  b(i) + u(i)*(-2*lambda*(1-theta) + 1);
    
    if i~=1 && i~=Mx-1
        b(i) = b(i) + u(i-1)*(lambda*(1-theta)) + u(i+1)*(lambda*(1-theta));
    elseif i==1
        b(i) = b(i) + g(1)*(lambda*(1-theta)) + u(i+1)*(lambda*(1-theta)) + g(1)*lambda*theta;
    elseif i==Mx-1
        b(i) = b(i) + u(i-1)*(lambda*(1-theta)) + g(Mx+1)*(lambda*(1-theta)) + g(Mx+1)*lambda*theta;
    end
end
b = -b;

for j = 1:T
    u = (A \ b);
    U(j,:) = u';
    %reinitialization
    b = zeros(Mx-1,1);

    for i = 1:Mx-1
        b(i) =  b(i) + u(i)*(-2*lambda*(1-theta) + 1);
        
        if i~=1 && i~=Mx-1
            b(i) = b(i) + u(i-1)*(lambda*(1-theta)) + u(i+1)*(lambda*(1-theta));
        elseif i==1
            b(i) = b(i) + g(1)*(lambda*(1-theta)) + u(i+1)*(lambda*(1-theta))+ g(1)*lambda*theta;
        elseif i==Mx-1
            b(i) = b(i) + u(i-1)*(lambda*(1-theta)) + g(Mx+1)*(lambda*(1-theta)) + g(Mx+1)*lambda*theta;
        end
    end
    b = -b;
    j/T*100
end
figure();
plot(g);
title('Initial condition g');
axis([0 Mx+1 0 1.2]);

figure();
for j=1:size(U,1)
    hold on;
    plot(U(j,:));
    axis([0 Mx+1 0 1.2]);
end