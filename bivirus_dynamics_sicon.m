%% This code simulates the deterministic bivirus model with higher order interactions

clear variables
close all
clc

tic

%Declare global variables for use in the simulation
global D1 D2 A1 A2 B1 B2 b1 b2 n

%% Parameter setup

n = 4;   %Number of nodes in the network

%By selecting the case number for the different parameters, you can quickly
%cycle between different scenarios.

recov = 1; %Diagonal recovery matrix case number
switch recov
    case 1
        D1 = 1.*eye(n);
        D2 = 1.*eye(n);
    case 2
        D1 = 1.*eye(n);
        D2 = 2.*eye(n);
end

infect = 4; %A^k captures the pairwise (1st order) infection effects for virus k.
switch infect
    case 1
        A1 = [2 0 0 0;1 1 1 0;0 1 1 0;0 0 0 0];
        A2 = [1 0 0 0;0 1 1 0;1 1 2 0;0 0 0 0];
        A1 = 2.*A1; A2 = 5.*A2;

    case 2  %Strange case to be studied later...multistability with copied viruses
        A1 = [1 0 0 0;1 1 1 0;0 1 1 0;0 0 1 0];
        A2 = [1 0 0 0;1 1 1 0;0 1 1 0;0 0 1 0];
        A1 = 2.*A1; A2 = 2.*A2;



    case 3  %Proposition 4.1
        A1 = [2 0 0 0;1 1 1 0;0 1 1 0;0 0 0 0];
        A2 = [1 0 0 0;0 1 1 0;1 1 2 0;0 0 0 0];
        A1 = 0.2.*A1; A2 = 0.2.*A2;

    case 4 % Proposition 4.3
        A1 = [2 0 0 0;1 1 1 0;0 1 1 0;0 0 0 0];
        A2 = [1 0 0 0;0 1 1 0;1 1 2 0;0 0 0 0];
        A1 = 2.*A1; A2 = 2.*A2;

    case 5 %Theorem 5.2
        A1 = [2 0 0 0;1 1 1 0;0 1 1 0;0 0 0 0];
        A2 = [1 0 0 0;0 1 1 0;1 1 2 0;0 0 0 0];
        A1 = 2.5.*A1; A2 = 2.*A2;

end


infect_hoi = 4;  
%B^k captures the HOI (2nd order) interaction terms for virus k. B1(:,:,i)
%captures the HOI represented by the matrix B^1_i in the paper, and
%similarly, B2(:,:,i) is the matrix B^2_i in the paper.

%The b_k correspond to the \beta_2^k parameters, viz the infection strength
%of the HOI interactions
switch infect_hoi
    case 1
        b1 = 0.2; b2 = 0.2;  

        B1 = zeros(n,n,n); B2 = B1;
        B1(:,:,1) = [1 0 0 0;0 0 1 0;0 1 0 1;0 0 1 0];
        B1(:,:,4) = [0 1 0 0;1 0 0 0;0 0 0 0;0 0 0 1];

        B2(:,:,1) = [1 0 0 0;0 0 1 1;0 1 0 0;0 1 0 0];
        B2(:,:,4) = [0 0 1 0;0 0 0 0;1 0 0 0;0 0 0 1];


    case 2 %Strange case to be studied later...multistability with copied viruses
        b1 = 5; b2 = 5;

        B1 = zeros(n,n,n); B2 = B1;
        B1(:,:,1) = [0 1 0 1;1 0 0 0;0 0 0 0;1 0 0 0];
        B1(:,:,4) = [0 0 0 1;0 0 0 1;0 0 0 0;1 1 0 0];

        B2(:,:,1) = [0 0 1 1;0 0 0 0;1 0 0 0;1 0 0 0];
        B2(:,:,4) = [0 0 0 1;0 0 0 0;0 0 0 1;1 0 1 0];

        B1(:,:,1) = [0 0 1 1;0 0 0 0;1 0 0 0;1 0 0 0];
        B1(:,:,4) = [0 0 0 1;0 0 0 0;0 0 0 1;1 0 1 0];

    case 3 %Proposition 4.1
        b1 = 10; b2 = 5;

        B1 = zeros(n,n,n); B2 = B1;
        B1(:,:,1) = [1 0 0 0;0 0 1 0;0 1 0 1;0 0 1 0];
        B1(:,:,4) = [0 1 0 0;1 0 0 0;0 0 0 0;0 0 0 1];

        B2(:,:,1) = [1 0 0 0;0 0 1 1;0 1 0 0;0 1 0 0];
        B2(:,:,4) = [0 0 1 0;0 0 0 0;1 0 0 0;0 0 0 1];

        bk = [1;0;0;1]; v = bk'*B1(:,:,1)*bk; u = A1*bk;

    case 4 %Proposition 4.3
        b1 = 0.5; b2 = 0.5;

        B1 = zeros(n,n,n); B2 = B1;
        B1(:,:,1) = [1 0 0 0;0 0 1 0;0 1 0 1;0 0 1 0];
        B1(:,:,4) = [0 1 0 0;1 0 0 0;0 0 0 0;0 0 0 1];

        B2(:,:,1) = [1 0 0 0;0 0 1 1;0 1 0 0;0 1 0 0];
        B2(:,:,4) = [0 0 1 0;0 0 0 0;1 0 0 0;0 0 0 1];

    case 5 %Theorem 5.2
        b1 = 2; b2 = 2.5;

        B1 = zeros(n,n,n); B2 = B1;
        B1(:,:,1) = [1 0 0 0;0 0 1 0;0 1 0 1;0 0 1 0];
        B1(:,:,4) = [0 1 0 0;1 0 0 0;0 0 0 0;0 0 0 1];

        B2(:,:,1) = [1 0 0 0;0 0 1 1;0 1 0 0;0 1 0 0];
        B2(:,:,4) = [0 0 1 0;0 0 0 0;1 0 0 0;0 0 0 1];

end

%This computes some quantities that are relevant to checking of conditions
%in the SIADS paper for the hypergraph to be strongly connected. That is,
%the following is relevant to Corollary 2.5
v = ones(1,n);
R1 = [v*B1(:,:,1);v*B1(:,:,2);v*B1(:,:,3);v*B1(:,:,4)];
R2 = [v*B2(:,:,1);v*B2(:,:,2);v*B2(:,:,3);v*B2(:,:,4)];
Z1 = A1+R1;
Z2 = A2+R2;
G1 = digraph(Z1); g1 = conncomp(G1);
G2 = digraph(Z2); g2 = conncomp(G2);


IC = 3; %Initial conditions case number. Allows for switching between different types of randomised IC generation.
switch IC
    case 1  %IC in the interior 
        b = rand(n,1);  a = rand(n,1); c = rand(n,1);
        x1 = a./(b+a+c); x2 = b./(a+b+c);
    case 2  %Close to one extreme edge of the boundary
        %         x1 = rand(n,1);  x2 = 0.00.*ones(n,1);
        x2 = 0.8.*rand(n,1);  x1 = 0.2.*ones(n,1);

    case 3  %Close to the origin
        x2 = 0.1.*rand(n,1); x1 = 0.1.*rand(n,1);
    case 4
        %         x1 = 0.95+0.49.*rand(n,1);  x2 = 0.0001.*ones(n,1);
        %                 x2 = 0.95+0.49.*rand(n,1);  x1 = 0.000000.*ones(n,1);
    case 5
        eta = 0.1;
        x1 = (0.5*eta).*ones(n,1); x2 = (1-eta).*ones(n,1);  %Bottom right
        %         x1 = (1-eta).*ones(n,1); x2 = (0.5*eta).*ones(n,1);  %Top left

end
x0 = [x1', x2']';

%% Simulate bivirus dynamics

s1 = max(real(eig(-D1+A1)));   disp(['s1 = ',num2str(s1)])
s2 = max(real(eig(-D2+A2)));   disp(['s2 = ',num2str(s2)])


tspan = [0 2000];  %Change 2nd term to set the simulation time
[t,x] = ode45(@bivirus_hoi,tspan,x0'); t = t'; x = x';   %Bivirus simulation

%Records the final values of virus 1 and virus 2, for convenient checking
%of the equilibria values.
x_final = x(:,end);
x1_final = x_final(1:n); x2_final = x_final(n+1:end);

%This is used for examining the stability properties at boundary equilibria
%only
p3 = max(eig(-D1+(eye(n)-diag(x2_final))*A1)); disp(['rho[(D^1)^{-1}(I-Z^2)B^1] = ',num2str(p3)])
p4 = max(eig(-D2+(eye(n)-diag(x1_final))*A2)); disp(['rho[(D^2)^{-2}(I-Z^1)B^2] = ',num2str(p4)])


%% Plot simulation output

figure
hold on
axis([0 tspan(2) 0 1])
xlab = xlabel('Time,  t');
ylab = ylabel('Fraction of Infected, x_i^k(t)');
p1 = plot(t,x(1:n,:),':b','LineWidth',1.5);
p2 = plot(t,x(n+1:end,:),'r','LineWidth',1.5);
leg = legend([p1(1) p2(1)],'Virus 1', 'Virus 2','box','off');
set(leg,'FontSize',12,'Location','SouthEast')
set(xlab,'FontSize',12)
set(ylab,'FontSize',12)

toc
