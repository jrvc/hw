%   Practical Excercise for the Galerkin approximation problem given by
%           a(u,phi) = (f,phi)   for all phi in H^1(0,2)
%   with a(u,phi) = (grad(u),grad(u)) + (u,phi)
clear all
clc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a) Calculate the STIFFNESS MATRIX A_h for the basis {1, x, ..., x^n}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 16;
%    the entries of A_h are a(phi_i,phi_j) where phi_k is the kth elem. of the basis
A_h = zeros(n+1);
for i = 0:n
    for j = i:n
        if (i == 0)
            A_h(i+1,j+1) = 2^(i+j+1) * ((i+j+1)^(-1));
        else
            A_h(i+1,j+1) = 2^(i+j-1) * (i*j*(i+j-1)^(-1)) + 2^(i+j+1) * ((i+j+1)^(-1));
        end
        A_h(j+1,i+1) = A_h(i+1,j+1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% NUMERICAL INTEGRATION SEEMS TO MAKE IT VERY MUCH MORE UNSTABLE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 0:n
%     for j = 0:n
%         fun1 = @(x) (i)*(j) * x.^(i+j-2);
%         int1 = integral(fun1,0,2);
%         fun2 = @(x) x.^((i+1)*(j+1));
%         int2 = integral(fun2,0,2);
%         A_h(i+1,j+1) = int1 + int2;
%         A_h(j+1,i+1) = A_h(i+1,j+1);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b) Calculate the l2-CONDITION NUMBER of the matrix A_h for for n = 1,...,16
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conds = zeros(n,1);
for k = 1:n
    conds(k) = cond(A_h(1:(k+1),1:(k+1)),2);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  c) CALCULATE ENTRIES OF b_h W.R.T. THE MONOMIAL BASIS for the right
%     hand side chosen s.t. the solution, u, is given by:
%           u(x) = cos (8*pi*x)
%     for n = 1,..,16
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(n+1,1);
bh= zeros(n+1,1);
for i=0:n
    fun3 = @(x) (-8) * (i) * x.^(i-1) * pi .* sin(8 * pi * x);
    fun4 = @(x) x.^(i+1) .* cos(8 * pi * x);
    b(i+1) = integral(fun3,0,2) + integral(fun4,0,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% PROFESSOR'S CODE FOR b_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r=2;
% b2 = zeros(n+1,1);
% syms variable;
% f = cos(8*pi*variable);
% for i=0:n
%     b2(i+1) = int( f .* variable^(i-1),variable,0,r);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d) Solve A_h * u_h = b_h
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Compute the coeficients 
    u_i = A_h \ b;
    X = [0:.001:2]; 
    % Compute the solution u_h as the combination of the elements of the
    % basis over the grid
    u_h = zeros(size(X));
    for j=1:1:n+1
        u_h = u_h + u_i(j) .* X.^(j-1);
    end

    

       
    % Compare the results with the real solution
    u_star = cos(8*pi .* X) ;
   
    plot(X, u_star)
    hold on
    plot(X, u_h)
    grid on
    title(' cos(8piX) and num. approx. polinomial basis')
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%  IGNORE THIS NEXT PART ~ NOT TRUE ANYMORE %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Can also be solved by Cholesky decomposition because the eigenvalues 
% are all positive

% if (all(eig(A_h) > 0))
%     %for n=1:17
%     L = chol(A_h); % produces L lower triangular s.t. L*L' = A
%     % First solve L*v = b for v = L'*u
%     v = L \ b;
%     v2 = L \ b2;
%     % Then solve L'*u = v
%     u = L' \ v;
%     u2 = L' \ v2;



