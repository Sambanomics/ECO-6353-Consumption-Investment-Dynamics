%%% Cooper, R., Haltiwanger, J., and Power, L. (1999)
%%% Machine Replacement and the Business Cycle: Lumps and Bumps
%%% (Replication Assignment)

%% Setting up The Value Function:

% V(k,A,e) = max[V^r(k,A,e), V^n(k,A,e)]
% where V^n(k,A,e) = Aef(k) + BEV(pk,A',e')
% V^r(k,A,e) = Aef(k)lamda - F + BEV(1,A',e')

% Parameter Value Allocation:
AH = 1.25;
AL = 0.75;
A = [AH,AL];
delta = 0.1;
beta = 0.9;
mu = 1;
rho = 1-delta;
F = 0.2;
lambda = 0.75;
N = 20;
R = 6; 
P=[0.9,0.1;0.1,0.9];

% Epsilon Grid:
e_min = 0.4;
e_max = 1.6;

e = linspace(e_min,e_max,N)';

% Capital Stock:

k = [1, cumprod(rho*ones(1,R - 1))];


%% 1. Initialzing the objects:

% (a) Creating a Three Dimensional Array:
s = e*k;
rev = s*A(1);
rev(:,:,2)= s*A(2);

% (b), (c) and (d) 3D guess for the Value Functions:

V_R = ones(size(rev));  
V_N = zeros(size(rev));  
V0 = max(V_R,V_N);       
 
tol = 10^(-9);           

%% 2. Writing the Value Function Iteration Loop: 
dif=1;
count=0;
maxits=10^9;

tic
while dif>tol && count<maxits
    z0 = zeros(size(rev));  
    for i=1:2
        V_R(:,:,i)= rev(:,:,i)*lambda -F + beta*(P(i,1)*repmat(V0(:,1,1),1,R)+P(i,2)*repmat(V0(:,1,2),1,R));
        for j=1:R-1
            V_N(:,j,i)= rev(:,j,i)+beta*(P(i,1)*V0(:,j+1,1)+P(i,2)*V0(:,j+1,2));
        end
    end
    V1 = max(V_R,V_N);
    z0(V1==V_R)=1;
    dif=max(max(max(abs(V1-V0))));
    count=count+1;
    V0=V1;
end
toc

%% 3. Examination of the Policy Function:

R_new = 7;

k1 = [1, cumprod(rho*ones(1,R_new - 1))];

% (a) Creating a Three Dimensional Array:           

s1 = e*k1;
rev1 = s1*A(1);
rev1(:,:,2)= s1*A(2);

% (b), (c) and (d) Three Dimensional Guess for the Value Functions:

V_R1 = ones(size(rev1));  
V_N1 = zeros(size(rev1));  
V00 = max(V_R1,V_N1);       
 
tol = 10^(-9);           


dif=1;
count=0;
maxits=10^9;

tic
while dif>tol && count<maxits
    z01 = zeros(size(rev1));  
    for i=1:2
        V_R1(:,:,i)= rev1(:,:,i)*lambda -F + beta*(P(i,1)*repmat(V00(:,1,1),1,R_new)+P(i,2)*repmat(V00(:,1,2),1,R_new));
        for j=1:R_new-1
            V_N1(:,j,i)= rev1(:,j,i)+beta*(P(i,1)*V00(:,j+1,1)+P(i,2)*V00(:,j+1,2));
        end
    end
    V11 = max(V_R1,V_N1);
    z01(V11==V_R1)=1;
    dif=max(max(max(abs(V11-V00))));
    count=count+1;
    V00=V11;
end
toc

%% 4. Plotting the Policy Function Using a "spy" Plot:

figure(1)
subplot(2,1,1)
spy(z0(:,:,1)')
xlabel('Firm Productivity')
ylabel('Time Since Last Replacement')
title('Replacement in Low Aggregate Productivity')

subplot(2,1,2)
spy(z0(:,:,2)')
xlabel('Firm Productivity')
ylabel('Time Since Last Replacement')
title('Replacement in High Aggregate Productvity')

% It can be inferred from the policy function that when the capital becomes obselete/ fully depreciated (reaches
% the cutoff), firms engage in investment for capital replacement and more inclined towards
% investing when the vintage of capital is older than 5 years. 
% Therefore, firms exhibits higher inclination to invest in high aggregate productivity
% state than low aggregagate productivity state in a scenario of high
% idiosycratic risk. (Drawing from the results in Cooper, R., Haltiwanger, J., and Power, L. (1999)

%% 5. Plotting the Hazard Function of Capital Replacement:

H = zeros(R,2);

for i=1:2
    for j=1:R
        H(j,i)=sum(z0(:,j,i))/N;
    end
end

time = 1:R ;

figure(2)
plot(time,H(:,1), 'r',time,H(:,2),'g')
title('Hazard Function for Machine Replacement')
xlabel('Time Since Last Replacement')
ylabel('Probability of Replacement')
legend('Low state','High State','Location','Best')
xlim([1 R])
ylim([0 1.05])

%% 6. Simuulating a Time Series and Plotting the Sample Paths:

ts=40;
T=160;
Time=(1:ts);
E=randi(N,1,ts);
esim=e(E);

AT=zeros(1,T);
AT(1)=2;





rand_num = rand;
if rand_num <= P(2)
    AT(i) = 2; 
else
    AT(i) = 1; 
end

AT = max(1, min(AT, numel(A)));

Asim=A(AT);

Y=zeros(1,ts);
K=zeros(1,ts+1);
K(1)=1;
sim_K=zeros(1,ts+1);
sim_K(1)=1;

for i=1:ts
    simoutput(i)=K(K(i))*Asim(i)*esim(i)*(1-z0(E(i),K(i),AT(i)))+z0(E(i),K(i),AT(i))*(lambda*K(K(i))*Asim(i)*esim(i)-F);
    if z0(E(i),K(i),AT(i))==1
        K(i+1)=1;
        sim_K(i+1)=K(K(i+1));
    else
        K(i+1)=K(i)+1;
        sim_K(i+1)=K(K(i+1));
    end
end

figure(3)
plot(Time,sim_K(1:ts),'r--',Time,simoutput,'g-')
title('Simulated Capital')
xlabel('Time')
ylabel('Evolution of the Firm')
ylim([0 max(simoutput)+0.1])
legend('Capital','Output','location','Best')


%% 7. Replicating Figure 3 in the Paper (Assuming that A is Fixed):

tt=50;
n=6;                     
ir=zeros(tt,1);
w=zeros(N,tt+1);
w(:,1)=ones(N,1);        
A_index=1;               

for i=1:tt
    for j=1:n
    ir(i)=ir(i)+H(j,A_index)*w(j,i);
    end
    w(1,i+1)=ir(i);
    for j=2:n
        w(j,i+1)=(1-H(j-1,A_index))*w(j-1,i);
    end
end

ir=ir/n;

figure(4)
plot(1:tt,ir,'red')
title('Convergence without Aggregate Shocks - Baseline Parameters')
xlabel('Time Period')
ylabel('Investment Rate')

%% 8. Following the Markov Process to Replicate Figure 4 in the Paper:

ir2=zeros(T-1,1);
w1=zeros(n,T);
w1(:,1)=ones(n,1); 

for i=1:T-1
    for j=1:n
    ir2(i)=ir2(i)+H(j, AT(i))*w1(j,i);
    end
    w1(1,i+1)=ir2(i);
    for j=2:n
        w1(j,i+1)=(1-H(j-1,AT(i)))*w1(j-1,i);
    end
end

ir2=ir2/n;

figure(5)
[y,line1,line2]=plotyy(1:T-1,ir2,1:T,Asim);
line1.Color = 'red'; 
line2.Color = 'green';
ylabel(y(1),'Investment Rate',Color="blue")
ylabel(y(2),'Aggregate State',Color="blue")
xlabel('period','FontWeight','bold')
line2.LineStyle=':';
line2.Marker='*';
title('Aggregate Investment Fluctuations - Baseline Simulation')

%% 9. Comment on the Implications (From Figure 3 and 4 From the Paper):

% Figures 4 and 5 depicts the dynamic patterns of aggregate investment levels 
% as they evolve across various stages of capital distribution, free of
% aggregate uncertainty. 

% Furthermore, the transition from low to high productivity states brings  
% about situations of a sharp increase in investments activities. 
% These fluctuations is depicted by Figures 4 and 5, indicating a situation
% where firms engage in lumpy investments to replace obselete capital vintages.