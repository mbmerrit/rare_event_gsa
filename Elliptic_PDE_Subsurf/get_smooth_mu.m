function mu_smooth = get_smooth_mu(mu, pde_data, theta)


load bg_param

p = pde_data.p;
e = pde_data.e;
t = pde_data.t;
g = pde_data.g;
    

c=theta;
q=pdeintrp(p,t, mu); 
alpha = 1;        
[K,F]=assempde(b,p,e,t,c,alpha,q);
mu_smooth = K\F;
        


return 


% testing stuff
u = mu_smooth;
figure(1);
pdeplot(p,[],t,'XYData',u,'Contour','on','XYStyle','flat',...
         'ColorBar','off','ColorMap','jet');

figure(2);
pdeplot(p,[],t,'XYData',u,'XYStyle','interp',...
         'ZData',u,'ZStyle','continuous',...
         'ColorBar','off');
colormap jet
