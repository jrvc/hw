% Solve the Dirichlet problem
%       u'' + b * u' + c * u = f
%                       u(0) = u(1) = 0
%

clear all
clc
close all

%rng(123)
%b = rand(1)
b = 22;
b_f=@(x) sin(x);
c = 105;
c_f=@(x)(0);
f = @(x)-exp(x);

T=1;

% set initial data
n = 100;
h = 1/n;
l = - 1/ h^2 - b/(2*h);
% grid
y = 0:h:1;
 
l_f=-1/h^2.*ones(1,n-1)-1/2/h.*feval(b_f,y(1:n-1));
d = 2/h^2 + c;
d_f=2/h^2.*ones(1,n-1)+feval(c_f,y(1:n-1));
r = - 1 / h^2 + b/(2*h);
r_f=-1/h^2.*ones(1,n-1)+1/2/h.*feval(b_f,y(1:n-1));
%boundary vconditions
g0 = 0;
g1 = 0;

%tridiagonal matrix
L_h = diag([1,d.*(ones(1,n-1)),1]) + diag([0, r.*ones(1,n-1)],1) + diag([l.*ones(1,n-1),0],-1);

L_f=diag([1,d_f,1],0)+diag([0,r_f],1)+diag([l_f,0],-1);


%f_h vector
f_h = [g0, feval(f,y(1:n-1)), g1];

% solve L_h u_h = f_h by different methods

% 1) Direct methods
if T==1
    
    u_h = L_h \ f_h';
    plot(y,u_h)
    grid on
    
else
    
    u_f = L_f \ f_h';
    figure()
    plot(y,u_f)
    grid on
end
% 2) Thomas algoritm
if T==1
    u_h2 = tridiag([1,d.*(ones(1,n-1)),1], [0,l.*ones(1,n-1),0], [0,r.*ones(1,n-1),0], f_h);
    figure()
    plot(y,u_h2)
    title('Thomas');
    grid on
else
    u_f2 = tridiag([1,d_f,1], [0,r_f,0], [0,l_f,0], f_h);
    figure()
    plot(y,u_f2)
    title('Thomas');
    grid on
end

cond(L_h)
% CG method
if T==1
    if (all(eig(L_h))>0)
        L_h = sparse(L_h);
        u_h3 = pcg(L_h,f_h',1e-7,2*n);
    else
        fprintf('Not possible to use CG')
    end
    figure()
    plot(y,u_h3)
    title('CG');
    grid on
else
    
    if (all(eig(L_f))>0)
        L_f = sparse(L_f);
        u_f3 = pcg(L_f,f_h',1e-7,1e+6);
    else
        fprintf('Not possible to use CG')
    end
    figure()
    plot(y,u_f3)
    title('CG');
    grid on
end
%GMRES method
if T==1
    
    u_h4 = gmres(L_h,f_h');
    figure()
    
    plot(y,u_h4)
    title('GMRES');
    grid on
    
else
    u_f4=gmres(L_f,f_h');
    figure()
    plot(y,u_f4)
    title('GMRES');
    grid on
end

function dydx = twoode(x,y)
dydx = [ y(2); -105*(y(1))-exp(x)-22*y(2) ];

function res = twobc(ya,yb)
res = [ ya(1); yb(1) ];

solinit = bvpinit(linspace(0,1),[0 0]);

sol = bvp4c(@twoode,@twobc,solinit);
figure()
x = linspace(0,4);
y = deval(sol,x);
plot(x,y(1,:))

    title('Shreck');
    grid on
