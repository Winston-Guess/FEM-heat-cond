clear all ; close all ; clc ;

% Set the number of elements
n_el = 6; n_np = 2*n_el + 1 ;

% Set the properties of the material and problem
A_e = 1 ; k_e = 0.45; 

%
s=500000000;
sx = inline('(500000000.*x)','x');

% Find the x-values for the nodal points 
x = linspace(0 , 0.02 , n_np ) ;

% Define the boundary conditions
T_A = 90; T_B = 205;

% Define the element stiffness matrix based on 
% variables defined in for loop
x_xi1 = -1/sqrt(3);
x_xi2 = 1/sqrt(3);
Ke_11 = inline('2*Ae*ke/le *((2 * x_xi1 -xe2 -xe3)^2 + (2 * x_xi2 -xe2-xe3)^2)'...
            ,'Ae','ke','le','x_xi1','x_xi2','xe2','xe3');
Ke_12 = inline('2*Ae*ke/le *(-2*(2 * x_xi1 -xe2 -xe3)*(2 * x_xi1 -xe1 -xe3) -2*(2 * x_xi2 -xe2 -xe3)*(2 * x_xi2 -xe1 -xe3))'...
            ,'Ae','ke','le','x_xi1','x_xi2','xe1','xe2','xe3');
Ke_13 = inline('2*Ae*ke/le *((2 * x_xi1 -xe2 -xe3)*(2 * x_xi1 -xe1 -xe2) + (2 * x_xi2 -xe2 -xe3)*(2 * x_xi2 -xe1 -xe2))','Ae'...
            ,'ke','le','x_xi1','x_xi2','xe1','xe2','xe3');
Ke_22 = inline('2*Ae*ke/le *(4*(2 * x_xi1 -xe1 -xe3)^2 + 4*(2 * x_xi2 -xe1 -xe3)^2)'...
            ,'Ae','ke','le','x_xi1','x_xi2','xe1','xe3');
Ke_23 = inline('2*Ae*ke/le *(-2*(2 * x_xi1 -xe1 -xe3)*(2 * x_xi1 -xe1 -xe2) -2*(2 * x_xi2 -xe1 -xe3)*(2 * x_xi2 -xe1 -xe2))'...
            ,'Ae','ke','le','x_xi1','x_xi2','xe1','xe2','xe3');
Ke_33 = inline('2*Ae*ke/le *((2 * x_xi1 -xe1 -xe2)^2 + (2 * x_xi2 -xe1 -xe2)^2)'...
            ,'Ae','ke','le','x_xi1','x_xi2','xe1','xe2');

% Define the calculation of element external body matrix
fe_11 = inline('1/le*( ((x_xi1-xe2)*(x_xi1-xe3))^2 + ((x_xi2-xe2)*(x_xi2-xe3))^2)'...
    ,'le','x_xi1','x_xi2','xe2','xe3');
fe_12 = inline('1/le*(-2*(x_xi1-xe2)*(x_xi1-xe3)^2*(x_xi1-xe1) -2*(x_xi2-xe2)*(x_xi2-xe3)^2*(x_xi2-xe1))'...
    ,'le','x_xi1','x_xi2','xe1','xe2','xe3');
fe_13 = inline('1/le*( (x_xi1-xe2)^2*(x_xi1-xe3)*(x_xi1-xe1) + (x_xi2-xe2)^2*(x_xi2-xe3)*(x_xi2-xe1))'...
    ,'le','x_xi1','x_xi2','xe1','xe2','xe3');
fe_22 = inline('1/le*( 4*((x_xi1-xe1)*(x_xi1-xe3))^2 + 4*((x_xi2-xe1)*(x_xi2-xe3))^2)'...
    ,'le','x_xi1','x_xi2','xe1','xe3');
fe_23 = inline('1/le*( -2*(x_xi1-xe1)^2*(x_xi1-xe2)*(x_xi1-xe3) -2*(x_xi2-xe1)^2*(x_xi2-xe2)*(x_xi2-xe3))'...
    ,'le','x_xi1','x_xi2','xe1','xe2','xe3');
fe_33 = inline('1/le*( ((x_xi1-xe1)*(x_xi1-xe2))^2 + ((x_xi2-xe1)*(x_xi2-xe2))^2 )'...
    ,'le','x_xi1','x_xi2','xe1','xe2');

% Define calculation of the gather index matrix
Le_indx = transpose([1 4:2:n_np-1; 3 5:2:n_np; 4 6:2:n_np-1 2]);

% Initialize the global stiffness matrix and external force matrix
K = zeros(n_np , n_np) ; f = zeros(n_np , 1) ;

% Initialize the matrix of element lengths
l_e = zeros(n_el,1); 

% Loop over elements
for e = 1:n_el         
    
    % Get the gather index matrix for current element in 'for' loop
    
    % Get the length of the current element in 'for' loop
    l_e(e) = x(2*e+1) - x(2*e-1);
    xe1 = x(2*e-1);
    xe2 = x(2*e);
    xe3 = x(2*e+1);
      
    % Use current gather index matrix to scatter and add Ke into K
    Ke11 = Ke_11(A_e,k_e,l_e(e),x_xi1,x_xi2,xe2,xe3);
    Ke12 = Ke_12(A_e,k_e,l_e(e),x_xi1,x_xi2,xe1,xe2,xe3);
    Ke13 = Ke_13(A_e,k_e,l_e(e),x_xi1,x_xi2,xe1,xe2,xe3);
    Ke22 = Ke_22(A_e,k_e,l_e(e),x_xi1,x_xi2,xe1,xe3);
    Ke23 = Ke_23(A_e,k_e,l_e(e),x_xi1,x_xi2,xe1,xe2,xe3);
    Ke33 = Ke_33(A_e,k_e,l_e(e),x_xi1,x_xi2,xe1,xe2);
    
    Ke = [Ke11 Ke12 Ke13;Ke12 Ke22 Ke23;Ke13 Ke23 Ke33];
    
    K(Le_indx(e,1:3),Le_indx(e,1:3)) = ...
                        K(Le_indx(e,1:3),Le_indx(e,1:3)) + Ke;
    
    s1=sx(xe1);    s2 = sx(xe2);    s3 = sx(xe3);
    fe11 = fe_11(l_e(e),x_xi1,x_xi2,xe2,xe3);
    fe12 = fe_12(l_e(e),x_xi1,x_xi2,xe1,xe2,xe3);
    fe13 = fe_13(l_e(e),x_xi1,x_xi2,xe1,xe2,xe3);
    fe22 = fe_22(l_e(e),x_xi1,x_xi2,xe1,xe3);
    fe23 = fe_23(l_e(e),x_xi1,x_xi2,xe1,xe2,xe3);
    fe33 = fe_33(l_e(e),x_xi1,x_xi2,xe1,xe2);
    
    fe = [fe11 fe12 fe13;fe12 fe22 fe23;fe13 fe23 fe33]*[s1 ; s2; s3];
%     % Use current gather index matrix to scatter and add fe into f
     f(Le_indx(e,1:3)) = f(Le_indx(e,1:3))+ fe;
end     

% Initialize the index matrix of unknown nodal tempretures
 dE_c_indx = 1 : 2;  dE_r_indx = 3 : n_np; 
 dF_indx = 3 : n_np ;

% Calculate unknown nodal tempretures using partitioning method
 dF = K(dF_indx, dF_indx) \ (f(dF_indx) -...
     (K(dE_r_indx, dE_c_indx)*[T_A ; T_B]));

% Construct the d matrix from known and calculated nodal tempretures
%d = [T_A ; T_B ; dF] ;

%Initialize the matrix of the gradient of tempreture along x
% grad_T = zeros( n_el );
% 
d = [T_A ; dF ;T_B]; 

grad_T = zeros(1,2);

Ne_1 = inline('2/(le^2)*(x-xe2)*(x-xe3)', 'le','x','xe2','xe3');
Ne_2 = inline('-4/(le^2)*(x-xe3)*(x-xe1)','le','x','xe1','xe3');
Ne_3 = inline('2/(le^2)*(x-xe1)*(x-xe2)', 'le','x','xe1','xe2');

Be_1 = inline('2/(le^2)*(2*x-xe2-xe3)', 'le','x','xe2','xe3');
Be_2 = inline('-4/(le^2)*(2*x-xe3-xe1)','le','x','xe1','xe3');
Be_3 = inline('2/(le^2)*(2*x-xe1-xe2)', 'le','x','xe1','xe2');

T_N=zeros(10*n_el,1);
xN=zeros(10*n_el,1);
xB=zeros(2,2);

figure

for e = 1 : n_el                
    % Get the current gradient of tempreture and insert into its matrix
    xe1 = x(2*e-1);
    xe2 = x(2*e);
    xe3 = x(2*e+1);
    
    c=1;
   
    for i = linspace (xe1,xe3,10);
        xN(10*(e-1)+c,1)=i;
        Ne1 = Ne_1(l_e(e),i,xe2,xe3);
        Ne2 = Ne_2(l_e(e),i,xe1,xe3);
        Ne3 = Ne_3(l_e(e),i,xe1,xe2);
        Ne = [Ne1 Ne2 Ne3];
        T_N(10*(e-1)+c,1) = Ne*d(2*e-1:2*e+1);
        c=c+1;
    end
    xB=[xe1;xe3];
    Be1_1 = Be_1(l_e(e),xe1,xe2,xe3);
    Be2_1 = Be_2(l_e(e),xe1,xe1,xe3);
    Be3_1 = Be_3(l_e(e),xe1,xe1,xe2);
    Be1 = [Be1_1 Be2_1 Be3_1];
    
    Be1_2 = Be_1(l_e(e),xe3,xe2,xe3);
    Be2_2 = Be_2(l_e(e),xe3,xe1,xe3);
    Be3_2 = Be_3(l_e(e),xe3,xe1,xe2);
    Be2 = [Be1_2 Be2_2 Be3_2];
    
    grad_T(1,1) = Be1 * d(2*e-1:2*e+1);
    grad_T(1,2) = Be2 * d(2*e-1:2*e+1);

    plot (xB, grad_T)
    hold on
end  

% Define plotting domain
x_plot = linspace(0, 0.02,40);

plot (x_plot , 115/.02 + s*(.02)^2/(6*A_e*k_e)-...
        ((s/(2*A_e*k_e)).*(x_plot.^2)),'Color','r')
title({'Heat Conduction Quadratic 1D FEM (10 Elements)';'vs Exact Solution'}...
            ,'FontWeight','normal','FontSize',12)
% Create xlabel
xlabel({'x (m)'},'FontWeight','normal','FontSize',12);
% Create ylabel
ylabel({'Tempreture gradient (m^{-1}K^{-1})'},'FontWeight','normal','FontSize',12);

% Plot
figure                                                                 
hold on ;
plot (xN , T_N);
plot (x_plot , 90 + (((s*(.02)^2)/(6*A_e*k_e))+...
    (115/0.02)).*x_plot - (s/(6*A_e*k_e)).* (x_plot.^3),'Color','r');
% legend ( 'FEM' , 'Exact' )                                          
title ( 'Quadratic 10 Elements vs Exact Solution',...
    'FontWeight','normal','FontSize',12 )
% Create xlabel
xlabel({'x (m)'},'FontWeight','normal','FontSize',12);
% Create ylabel
ylabel({'Tempreture (^{\circ}C)'},'FontWeight','normal','FontSize',12);
