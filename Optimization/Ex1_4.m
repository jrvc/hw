%
% Script for solving the numerical excercise in the 1st Optimization sheet
%
% MINIMAL TRIANGULATED GRAPHS
% Compute the minimal triangular graph over Omega := (0,1)x(0,1).
%
% Find q : Ömega -> R, s.t.  
%         a) q is picewise linear and continuos in a linear triangulation of Omega
%         b) q has the minimal surface graph among all picewise functionals
%         c) q satisfies the boundary conditions:
%               q(x1,x2) = 1/2 - |x2 - 1/2| for (x1,x2) in the boundary  
%
clear all
clc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subdivide Omega into triangles with corners Q1(x1', x2''),..., Qm(x1', x2'')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 10; % number of triangles is 2*n^2
h = 1/n; % size of the grid
% grid
[x,y] = meshgrid(0:h:1,0:h:1);
%triangulation
tri = delaunay(x,y);
tri  = sort(tri,2);
[Y,ind] = sort(tri(:,1));
B = tri(ind,:);
tri = B;
%prove if it works
z = zeros(n+1);
trisurf(tri,x,y,z)

% we need to be able to examine the coordinates of the vertices of the
% triangles, so we define the triangular grid as follows
P = [x(:,1),y(:,1)];
for i=2:n+1
    P = [P; x(:,i),y(:,i) ];
end


TR = triangulation(tri,P)
triplot(TR)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the functions q continuous on Omega and linear in each triangle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the values on the boundary are fixed, i.e., z(1,:),z(n+1,:),z(:,1),z(:,n+1)
z(:,1) = 1/2 - abs(y(:,1) - 1/2)
z(:,n+1) = 1/2 - abs(y(:,n+1) - 1/2)
trisurf(tri,x,y,z)

% we fix the initial conditions to be as random between 0 and 0.5
rng(12345);
z(:,2:n) = unifrnd(0,ones((n+1),(n-1)).*0.5);
trisurf(tri,x,y,z)

% Calculate the initial Areas
% Recall: the magnitud of the cross product between two vectors is equal to
%         the area that has the parallelogram that is described by its sum.

areas = zeros(1,(2*n^2));

% the instruction TR.Points(TR.ConnectivityList(k,:),:)  lets you examine 
% the coordinates of the vertices of the k-th triangle
for k = 1:(2*n^2)
    vertices = TR.ConnectivityList(k,:);
    triang = [TR.Points(vertices,:), [z(vertices)]' ];
    areas(k) = 0.5 * norm(cross(triang(2,:) - triang(1,:), triang(3,:)- triang(1,:)));
end
TotArea = sum(areas)


