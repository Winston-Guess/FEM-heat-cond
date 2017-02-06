clear all ; close all ; clc ;
npoints = 100;
% Set the number of elements
n_el = 18; n_np = 16 ;

% Set the properties of the material and problem
k_e = 5;
s = 6;

dE = [0;0;0;0;0;0;0;0];

ICA = [1 2 13            % Interconnectivity array
       2 3 13
       3 14 13
       3 4 14
       14 4 15
       4 16 15
       15 16 5
       15 5 6
       13 15 6
       13 14 15
       12 13 6
       12 6 7
       9 12 7
       9 7 8
       10 12 9
       10 11 12
       11 13 12
       11 1 13];
   
GA = [ 3  3.1340  3.5   4  4  3  1.5  0  0  0  1.5  1.5  3  3.5  3.5  4    % x
       0   0.5   .8660  1  2  2   2   2  1  0   0    1   1   1   1.5  1.5 ];      % y


A_e=inline('.5*((xe2*ye3-xe3*ye2)-(xe1*ye3-xe3*ye1)+(xe1*ye2-xe2*ye1))'...
    ,'xe1','xe2','xe3','ye1','ye2','ye3');

Be_11=inline('(ye2-ye3)','ye2','ye3');
Be_12=inline('(ye3-ye1)','ye1','ye3');
Be_13=inline('(ye1-ye2)','ye1','ye2');
Be_21=inline('(xe3-xe2)','xe2','xe3');
Be_22=inline('(xe1-xe3)','xe1','xe3');
Be_23=inline('(xe2-xe1)','xe1','xe2');

elx = zeros(n_el,4);                         % element x-coords var
ely = zeros(n_el,4);                         % element y-coords var

K = zeros(n_np , n_np);
f_Omega = zeros(n_np ,1);
f_Gamma = zeros(n_np ,1);

for e = 1:n_el

    for i=1:3
        elx(e,i) = GA(1,ICA(e,i));
        ely(e,i) = GA(2,ICA(e,i));
        
        if i==1
            elx(e,4) = GA(1,ICA(e,i));
            ely(e,4) = GA(2,ICA(e,i));
        end
    end

    xe1 = elx(e,1); ye1 = ely(e,1);
    xe2 = elx(e,2); ye2 = ely(e,2);
    xe3 = elx(e,3); ye3 = ely(e,3);

    Ae = A_e(xe1,xe2,xe3,ye1,ye2,ye3);

    Be11 = Be_11(ye2,ye3);
    Be12 = Be_12(ye1,ye3);
    Be13 = Be_13(ye1,ye2);
    Be21 = Be_21(xe2,xe3);
    Be22 = Be_22(xe1,xe3);
    Be23 = Be_23(xe1,xe2);

    Be = [Be11 Be12 Be13; Be21 Be22 Be23];

    Ke = k_e/4/Ae*transpose(Be)*Be;
        
    K(ICA(e,1:3),ICA(e,1:3)) = K(ICA(e,1:3),ICA(e,1:3)) + Ke;

    fe_Omega = s*Ae/3*[1;1;1];

    f_Omega(ICA(e,1:3)) = f_Omega(ICA(e,1:3)) + fe_Omega;

    hold on

end

f_Gamma(8:10) = 20*[.5;1;1];

f = f_Omega + f_Gamma;

dE_c_indx = 1 : 8;  dE_r_indx = 9 : n_np; 
dF_indx = 9 : n_np ;

% Calculate unknown nodal tempretures using partitioning method
dF = K(dF_indx, dF_indx) \ (f(dF_indx) - (K(dE_r_indx, dE_c_indx)*dE));

d = [dE ; dF];

for e = 1:n_el

patch(elx(e,1:3),ely(e,1:3),d(ICA(e,1:3)))
    
    %plot (elx(e,1:4), ely(e,1:4),'Color','w')
end
title('Heat Conduction Triangular Element FEM'...
            ,'FontWeight','normal','FontSize',12)
% Create xlabel
xlabel({'x (m)'},'FontWeight','normal','FontSize',12);
% Create ylabel
ylabel({'y (m)'},'FontWeight','normal'...
    ,'FontSize',12);
colorbar
legend
