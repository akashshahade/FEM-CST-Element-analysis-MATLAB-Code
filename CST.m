
%-------------------------------------------------------------------%
% Software by -  Akash S. Shahade
%
% Title - Program for FE Analysis of Plate using CST Element.
%-------------------------------------------------------------------%
clc;clear;
fprintf(' \nFEM ANALYSIS OF PLATE USING CST Element. \n\n');

%------------------------------------------------------%
% STEP 01 - PRE-PROCESSING
%------------------------------------------------------%

u = 0.25;          % Poisson's Ratio
E = 2*10^(11);     % Young's Modulus (N/m^2)
t = 0.015;         % Plate Thickness (m)

ele_nod=[1 3 4;1 2 3];                       % Elements are connected with these nodes
nod_coor=[0 0;0.750 0;0.750 0.500;0 0.500];  % Node coordinates
num_ele=size(ele_nod,1);                     % Number of elements
ele_dof=[1 2 5 6 7 8;1 2 3 4 5 6];           % D.O.F  associated with Nodes
num_nod=4;                                   % Number of Nodes
dof = 2;                                     % D.O.F per node

displacement = zeros(dof*num_nod,1);         % Zero Matrix for Displacement
force = zeros(dof*num_nod,1);                % Zero Matrix for Force
stiffness = zeros(dof*num_nod);              % Zero Matrix for Stiffness

s = E/(1 - (u^2)) ;
D = s*[1 u 0; u 1 0; 0 0 (1 - u)*0.5];       % D matrix 

J1 = [1 nod_coor(1,1) nod_coor(1,2);...      % Jacobian of Element 01
    1 nod_coor(3,1) nod_coor(3,2);...
    1 nod_coor(4,1) nod_coor(4,2)];

J2 = [1 nod_coor(1,1) nod_coor(1,2);...      % Jacobian of Element 02
    1 nod_coor(2,1) nod_coor(2,2);...
    1 nod_coor(3,1) nod_coor(3,2)];

A1 = 0.5 * det(J1);                          % Area of Element 01
A2 = 0.5 * det(J2);                          % Area of Element 02

% Element 01 [B] Matrix

B1 = [(nod_coor(3,2)-nod_coor(4,2)) 0 (nod_coor(4,2)-nod_coor(1,2)) 0 (nod_coor(1,2)-nod_coor(3,2)) 0;...
    0 (nod_coor(4,1)-nod_coor(3,1)) 0 (nod_coor(1,1)-nod_coor(4,1)) 0 (nod_coor(3,1)-nod_coor(1,1));...
    (nod_coor(4,1)-nod_coor(3,1)) (nod_coor(3,2)-nod_coor(4,2)) (nod_coor(1,1)-nod_coor(4,1))...
    (nod_coor(4,2)-nod_coor(1,2)) (nod_coor(3,1)-nod_coor(1,1)) (nod_coor(1,2)-nod_coor(3,2))];

% Element 02 [B] Matrix

B2 = [(nod_coor(2,2)-nod_coor(3,2)) 0 (nod_coor(3,2)-nod_coor(1,2)) 0 (nod_coor(1,2)-nod_coor(2,2)) 0;...
    0 (nod_coor(3,1)-nod_coor(2,1)) 0 (nod_coor(1,1)-nod_coor(3,1)) 0 (nod_coor(2,1)-nod_coor(1,1));...
    (nod_coor(3,1)-nod_coor(2,1)) (nod_coor(2,2)-nod_coor(3,2)) (nod_coor(1,1)-nod_coor(3,1))...
    (nod_coor(3,2)-nod_coor(1,2)) (nod_coor(2,1)-nod_coor(1,1)) (nod_coor(1,2)-nod_coor(2,2))];

%------------------------------------------------------%
% Stiffness matrix calculation & ASSEMBLY
%------------------------------------------------------%

for e=1:num_ele                               % For 1 to Number of elements
   
    if e == 1
        k = (B1.'*D*B1)*t*A1 ;
    else
        k = (B2.'*D*B2)*t*A2 ;
    end

   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);

    for i=1:6
        for j=1:6
                                              % Assembly of Global Matrix
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);

        end
    end
end

force = [0;0;0;0;50000;0;0;0];                % Force Vector


fprintf('Global Stiffness Matrix: \n');
disp(stiffness);
fprintf('\n Global Load Vector: \n');
disp(force);
fprintf('\n------------------------------------------------\n');

%------------------------------------------------------%
% Boundary Conditions
%------------------------------------------------------%

fixed_dof = [1 2 4 7 8];             % Fixed D.O.F.
k=stiffness;
k(fixed_dof,:)=[];                   % Eliminating Rows
k(:,fixed_dof)=[];                   % Eliminating Columns

f=force;
f(fixed_dof,:)=[];                   % Eliminating Rows


%------------------------------------------------------%
% STEP 02 - SOLVE
%------------------------------------------------------%

q = k\f ;

%------------------------------------------------------%
% STEP 03 - POST-PROCESSING
%------------------------------------------------------%

displacement = [0;0;q(1);0;q(2);q(3);0;0];   % Displacement Vector

reaction = stiffness*displacement - force;   % Calculate Reaction Forces

ele1_disp = [1 2 5 6 7 8];                   % Element 01 Displacements
q1 = displacement(ele1_disp,:);  

ele1_disp = [1 2 3 4 5 6];                   % Element 02 Displacements
q2 = displacement(ele1_disp,:);

Stress1 = D*B1*q1;                           % Stress in Element 01
Stress2 = D*B2*q2;                           % Stress in Element 02


q={'q1x';'q1y';'q2x';'q2y';'q3x';'q3y';'q4x';'q4y'};
Displacement_mm = displacement*1000;
T = table(q,Displacement_mm);
disp(T);                                     % Print Displacement

fprintf('\n------------------------------------------------\n');

R = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8'};
P = table(R,reaction);
disp(P);                                     % Print Reaction

fprintf('\n------------------------------------------------\n');

Element = [1;1;1;2;2;2];
stress = [Stress1;Stress2];
Stress_MPa = (stress)/10^6;
B = table(Element,Stress_MPa);
disp(B);                                     % Print Stresses

fprintf('END OF PROGRAM.\n');
