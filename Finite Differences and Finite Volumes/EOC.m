close all;

m = 20; %%% 5
for j=1:100/m
    v_direct_5(j) = u_h_100(m*j) - u_h_5(j);      
    v_thomas_5(j) = u_h2_100(m*j) - u_h2_5(j);    
    v_cg_5(j) = u_h3_100(m*j) - u_h3_5(j);        
    v_gmres_5(j) = u_h4_100(m*j) - u_h4_5(j);    
end

m = 10; %%% 10
for j=1:100/m
    v_direct_10(j) = u_h_100(m*j) - u_h_10(j);      
    v_thomas_10(j) = u_h2_100(m*j) - u_h2_10(j);    
    v_cg_10(j) = u_h3_100(m*j) - u_h3_10(j);        
    v_gmres_10(j) = u_h4_100(m*j) - u_h4_10(j);    
end

m = 5; %%% 20
for j=1:100/m
    v_direct_20(j) = u_h_100(m*j) - u_h_20(j);      
    v_thomas_20(j) = u_h2_100(m*j) - u_h2_20(j);    
    v_cg_20(j) = u_h3_100(m*j) - u_h3_20(j);        
    v_gmres_20(j) = u_h4_100(m*j) - u_h4_20(j);    
end

m = 4; %%% 25
for j=1:100/m
    v_direct_25(j) = u_h_100(m*j) - u_h_25(j);      
    v_thomas_25(j) = u_h2_100(m*j) - u_h2_25(j);    
    v_cg_25(j) = u_h3_100(m*j) - u_h3_25(j);        
    v_gmres_25(j) = u_h4_100(m*j) - u_h4_25(j);    
end


m = 2; %%% 50
for j=1:100/m
    v_direct_50(j) = u_h_100(m*j) - u_h_50(j);      
    v_thomas_50(j) = u_h2_100(m*j) - u_h2_50(j);    
    v_cg_50(j) = u_h3_100(m*j) - u_h3_50(j);        
    v_gmres_50(j) = u_h4_100(m*j) - u_h4_50(j);    
end

E_direct(1) = norm(v_direct_5);
E_direct(2) = norm(v_direct_10);
E_direct(3) = norm(v_direct_20);
E_direct(4) = norm(v_direct_25);
E_direct(5) = norm(v_direct_50);

E_thomas(1) = norm(v_thomas_5);
E_thomas(2) = norm(v_thomas_10);
E_thomas(3) = norm(v_thomas_20);
E_thomas(4) = norm(v_thomas_25);
E_thomas(5) = norm(v_thomas_50);

E_cg(1) = norm(v_cg_5);
E_cg(2) = norm(v_cg_10);
E_cg(3) = norm(v_cg_20);
E_cg(4) = norm(v_cg_25);
E_cg(5) = norm(v_cg_50);

E_gmres(1) = norm(v_gmres_5);
E_gmres(2) = norm(v_gmres_10);
E_gmres(3) = norm(v_gmres_20);
E_gmres(4) = norm(v_gmres_25);
E_gmres(5) = norm(v_gmres_50);

H = [1/5 1/10 1/20 1/25 1/50];

E_direct(1) = [];
E_thomas(1) = [];
E_cg(1) = [];
E_gmres(1) = [];

H(1) = [];

figure();
plot(log(H),log(E_direct),'o');
title('Direct');
grid;
X = [ones(length(log(H)'),1) log(H)'];
linear_reg_direct = X\log(E_direct)'
h = lsline
set(h(1),'color','r')

figure();
plot(log(H),log(E_thomas),'o');
title('Thomas');
grid;
p_thomas = log(E_thomas) / log(H)
X = [ones(length(log(H)'),1) log(H)'];
linear_reg_thomas = X\log(E_thomas)'
h = lsline
set(h(1),'color','r')

figure();
plot(log(H),log(E_cg),'o');
title('CG');
grid;
p_cg = log(E_cg) / log(H)
X = [ones(length(log(H)'),1) log(H)'];
linear_reg_cg = X\log(E_cg)'
h = lsline
set(h(1),'color','r')

figure();
plot(log(H),log(E_gmres),'o');
title('GMRES');
grid;
p_gmres = log(E_gmres) / log(H)
X = [ones(length(log(H)'),1) log(H)'];
linear_reg_gmres = X\log(E_gmres)'
h = lsline
set(h(1),'color','r')