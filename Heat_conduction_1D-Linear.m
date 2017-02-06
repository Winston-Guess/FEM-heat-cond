clear all ; close all ; clc ;

% Set the number of elements
n_el = 6; n_np = n_el + 1 ;

% Set the properties of the material and problem
A_e = 1 ; k_e = 0.45 ; 

%
s=500000000;
sx = inline('(500000000.*x)','x');

% Find the x-values for the nodal points 
x = linspace(0 , 0.02 , n_np ) ;

% Define the boundary conditions
T_A = 90; T_B = 205;

% Define the element stiffness matrix based on 
% variables defined in for loop
Ke = inline('(Ae * ke / le) .* [1 -1;-1 1]','Ae','ke','le');

% Define the calculation of element external body matrix
fe_Omega = inline('(le/6).*[2 1 ; 1 2]*[s1 ; s2] ','le','s1','s2') ;

% Define the calculation of the gradient of the shape function matrix
Be = inline('( 1/le ) .* [-1 1]' , 'le');

% Define calculation of the gather index matrix
Le_indx = transpose([1 3:n_np-1 n_np; 3 4:n_np 2]);

% Initialize the global stiffness matrix and external force matrix
K = zeros(n_np , n_np) ; f = zeros(n_np , 1) ;

% Initialize the matrix of element lengths
l_e = zeros(n_el,1); 

% Loop over elements
for e = 1:n_el         
    
    % Get the length of the current element in 'for' loop
    l_e(e) = x(e+1) - x(e);
    
    % Use current gather index matrix to scatter and add Ke into K
    K(Le_indx(e,1:2),Le_indx(e,1:2)) = ...
        K(Le_indx(e,1:2),Le_indx(e,1:2)) + Ke(A_e , k_e , l_e(e));
    
    x1=x(e);    x2=x(e+1);
    s1=sx(x1);    s2 = sx(x2);
    % Use current gather index matrix to scatter and add fe into f
    f(Le_indx(e,1:2)) = f(Le_indx(e,1:2))+ fe_Omega(l_e(e) , s1 , s2);
end     

% Initialize the index matrix of unknown nodal tempretures
dE_c_indx = 1 : 2;  dE_r_indx = 3 : n_np; 
dF_indx = 3 : n_np ;

% Calculate unknown nodal tempretures using partitioning method
dF = K(dF_indx, dF_indx) \ (f(dF_indx) - ...
                                (K(dE_r_indx, dE_c_indx)*[T_A ; T_B]));

% Construct the d matrix from known and calculated nodal tempretures
d = [T_A ; T_B ; dF] ;

%Initialize the matrix of the gradient of tempreture along x
grad_T = zeros(n_el,1);

figure
for e = 1 : n_el                       
    % Get the current gradient of tempreture and insert into its matrix
    grad_T(e) = Be(l_e(e)) * d(Le_indx(e,1:2));
    
    plot ( x(e:e+1),[grad_T(e) grad_T(e)])
    hold on ;
    % legend ( 'FEM' , 'Exact' )

end

d = [T_A ; dF ; T_B];

xi_e = 1/sqrt(3)*[-1 1];

% Define the calculation of the shape function matrix for error analysis
N_e = inline('( 1/le ) .*[xe_2-x x-xe_1]*[d1; d2]', ...
                                'le','x','xe_1','xe_2','d1','d2');
                            
U_ex = inline('90+(((s*(.02)^2)/(6*A_e*k_e))+(115/0.02))*x-(s/(6*A_e*k_e))*(x^3)',...
                'A_e','k_e','s','x');
eL2_2 = 0;
eL2_e_bar2=0;

%Loop for calculationg Error
for e = 1 : n_el
    eL2_e=0;
    eL2_e_bar=0;
    
    d1 = d(e);
    d2 = d(e+1);
    xe_1 = x(e);
    xe_2 = x(e+1);
    for x_i = xi_e
        xi = .5*(xe_1+xe_2)+.5*x_i*(xe_2-xe_1);
        Ned_xi = N_e (l_e(e),xi,xe_1,xe_2,d1,d2);
        Uex_xi = U_ex (A_e,k_e,s,xi);
        eL2_e=eL2_e+(Uex_xi - Ned_xi)^2;
        eL2_e_bar = eL2_e_bar+(Uex_xi - Ned_xi)^2/Uex_xi^2;
    end
    eL2_2= eL2_2 + eL2_e;
    eL2_e_bar2= eL2_e_bar2 + eL2_e_bar;
end

% Error in the L2 norm 
eL2=eL2_2^(1/2);
% Normalized Error in the L2 norm
eL2_bar=eL2_e_bar2^(1/2);

% Define plotting domain
x_plot = linspace(0, 0.02,40);

plot (x_plot , 115/.02 + s*(.02)^2/(6*A_e*k_e)-...
    ((s/(2*A_e*k_e)).*(x_plot.^2)),'Color','r')

title({'Heat Conduction Linear 1D FEM (10 Elements)';'vs Exact Solution'}...
            ,'FontWeight','normal','FontSize',12)
% Create xlabel
xlabel({'x (m)'},'FontWeight','normal','FontSize',12);
% Create ylabel
ylabel({'Tempreture gradient (m^{-1}K^{-1})'},'FontWeight','normal'...
    ,'FontSize',12);

% Plot
figure                                                                 
hold on ;                                                               
plot (x_plot , 90 + (((s*(.02)^2)/(6*A_e*k_e))+(115/0.02)).*x_plot -...
    (s/(6*A_e*k_e)).* (x_plot.^3),'Color','r')

plot( x , d)
% legend ( 'FEM' , 'Exact' )                                          
title ( 'Tempreture in 1D Heat Conduction as a function of length',...
    'FontWeight','normal','FontSize',12)
% Create xlabel
xlabel({'Length (m)'},'FontWeight','normal','FontSize',12);
% Create ylabel
ylabel({'Tempreture (*K)'},'FontWeight','normal','FontSize',12);
