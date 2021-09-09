clc;
clear;
close all;
format short g
%% Problem Definition


% Kpss= [0.1 50];
% T1=[0.2  1.5];
% T2=[0.02  0.15];
% T3=[0.2   1.5];
% T4=[0.02  0.15];

CostFunction=@(p)costsimb(p);        % Cost Function

nVar=3;       % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix
%       T1   T2  KPSS
VarMin=[0.2 0.02 0.1];         % Lower Bound of Variables
VarMax=[1.5 0.15 50];         % Upper Bound of Variables


%% PSO Parameters

MaxIt=100;      % Maximum Number of Iterations

nPop=20;        % Population Size (Swarm Size)
wmin=0.4;
wmax=0.9;
c1i=2.5;
c2i=0.2;
c1f=0.2;
c2f=2.5;
alfa =1;
R=0.25;
s=1;
Q=-1*(1/MaxIt)*log(exp(s));

for it=1:MaxIt
c1(it)=((c1f-c1i).*(it/MaxIt))+c1i;
c2(it)=((c2f-c2i).*(it/MaxIt))+c2i;
w(it)=wmax-(wmax-wmin).*(1+tanh(Q*(it-(MaxIt/2))))/2;
w(it)=wmin+rand.*w(it);
end


% Velocity Limits
VelMax=R*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Sol=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.Best.Sol=[];
particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    % Initialize Position
    particle(i).Position=random('unif',VarMin,VarMax);
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    [particle(i).Cost particle(i).Sol]=CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.Sol=particle(i).Sol;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

nfe=zeros(MaxIt,1);


%% PSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w(it)*particle(i).Velocity ...
            +c1(it)*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2(it)*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
   
        particle(i).Position2 = particle(i).Position + particle(i).Velocity;
        
        if particle(i).Position2 <particle(i).Position
            particle(i).Position = particle(i).Position -alfa*abs (particle(i).Velocity);
        elseif  particle(i).Position2 >particle(i).Position
           particle(i).Position  = particle(i).Position +alfa* abs(particle(i).Velocity); 
        else
             
        particle(i).Position = particle(i).Position2; 
        end
            
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        [particle(i).Cost particle(i).Sol] = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).Best.Sol=particle(i).Sol;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    % Store the best cost value
    BestCost(it)=GlobalBest.Cost;
    
    % Display Iteration information
    disp(['Iteration ' num2str(it)  ', Best Cost = ' num2str(BestCost(it))]);
    
    %w=w*wdamp;
    
end

%% Results

figure(1);
plot(BestCost,'LineWidth',2);
xlabel('Itiration');
ylabel('Best Cost');
num2str(GlobalBest.Position)

%%
Re=0;
Xe=0.5;
Vt=cos(15*pi/180)+sqrt(-1)*sin(15*pi/180);
Vinf=1.05;
H=3.25;
Tpdo=9.6;
KA=400;
TA=0.2;
Rs=0;
Xq=2.1;
Xd=2.5;
Xpd=0.39;
D=0;
ws=377;
T=3;
Time = 20;
dVref = 0.0;
dTm= 0.2;
tsim= 1.0005;     % 1.0005 <tsim< 20  sec
% ========================>>> Find Initial Value <<<========================
i=sqrt(-1);
IG=(Vt-Vinf)/(Re+sqrt(-1)*Xe);
Vt0=Vt;
gama=angle(IG);
teta=angle(Vt0);
delta0=angle(Vt0+(Rs+sqrt(-1)*Xq)*IG);
Idq= abs(IG)*(cos(pi/2+gama-delta0)+sqrt(-1)*sin(pi/2+gama-delta0));
Vdq=abs(Vt0)*(cos(pi/2+teta-delta0)+sqrt(-1)*sin(pi/2+teta-delta0));
Id0=real(Idq);
Iq0=imag(Idq);
Vd0=real(Vdq);
Vq0=imag(Vdq);
Epd0=Vd0+Rs*Id0-Xpd*Iq0;
Epq0=Vq0+Rs*Iq0+Xpd*Id0;
Efd0=Epq0+(Xd-Xpd)*Id0;
Vref=abs(Vt0)+Efd0/KA;
Omega=ws;
TM=Epq0*Iq0+(Xq-Xpd)*Id0*Iq0;

% ===========>>> Calculation K1… K6 ,values of  Hiffron-Flips Model <<<===========
DELTA=Re^2+(Xe+Xq)*(Xe+Xpd);
K1=(-1/DELTA)*(Iq0*Vinf*(Xpd-Xq)*((Xq+Xe)*sin(delta0)-Re*cos(delta0))+...
Vinf*((Xpd-Xq)*Id0-Epq0)*((Xpd+Xe)*cos(delta0)+Re*sin(delta0)));
K2=(1/DELTA)*(Iq0*DELTA-Iq0*(Xpd-Xq)*(Xq+Xe)-Re*(Xpd-Xq)*Id0+Re*Epq0);
K3=1/(1+((Xd-Xpd)*(Xq+Xe)/DELTA));
K4=(Vinf*(Xd-Xpd)/DELTA)*((Xq+Xe)*sin(delta0)-Re*cos(delta0));
K5=(1/DELTA)*((Vd0/abs(Vt0))*Xq*(Re*Vinf*sin(delta0)+Vinf*cos(delta0)*...
(Xpd+Xe))+(Vq0/abs(Vt0))*(Xpd*(Re*Vinf*cos(delta0)-Vinf*(Xq+Xe)*...
sin(delta0))));
K6=(1/DELTA)*((Vd0/abs(Vt0))*Xq*Re-(Vq0/abs(Vt0))*Xpd*(Xq+Xe))+...
(Vq0/abs(Vt0));
T1=GlobalBest.Position(1);
T2=GlobalBest.Position(2);
Kpss=GlobalBest.Position(3);

%% PLOT
% PSS= 1; %#ok Simulation with PSS
% simopt = simset('solver','ode45','SrcWorkspace','Current','DstWorkspace','Current');  % Initialize sim options
% [tout xout yout]=sim('SMIB_sim',[0 Time],simopt);  %#ok
% figure (4)
% subplot(3,1,1)
% plot(ScopeData2.time,ScopeData2.signals(1, 1).values,'LineWidth',1.5);hold on
% title('SMIB ')
% ylabel('\Deltas')
% subplot(3,1,2)
% plot(ScopeData2.time,ScopeData2.signals(1, 2).values,'LineWidth',1.5);hold on
% ylabel('{\Delta}{\delta}')
% subplot(3,1,3)
% plot(ScopeData2.time,ScopeData2.signals(1, 3).values,'LineWidth',1.5);hold on
% ylabel('{E_q}')
% xlabel('Time')
% hold on
% PSS= 0; % Simulation with PSS
% simopt = simset('solver','ode45','SrcWorkspace','Current','DstWorkspace','Current');  % Initialize sim options
% [tout xout yout]=sim('SMIB_sim',[0 Time],simopt);  %#ok
% figure (4)
% subplot(3,1,1)
% plot(ScopeData2.time,ScopeData2.signals(1, 1).values,'r','LineWidth',1.5)
% ylabel('\Deltas')
% subplot(3,1,2)
% plot(ScopeData2.time,ScopeData2.signals(1, 2).values,'r','LineWidth',1.5);legend('with PSS', 'without PSS')
% ylabel('{\Delta}{\delta}')
% subplot(3,1,3)
% plot(ScopeData2.time,ScopeData2.signals(1, 3).values,'r','LineWidth',1.5)
% ylabel('{E_q}')
% xlabel('Time')

PSS= 0; %#ok Simulation without PSS
simopt = simset('solver','ode45','SrcWorkspace','Current','DstWorkspace','Current');  % Initialize sim options
[tout xout yout]=sim('SMIB_sim',[0 Time],simopt);  %#ok
figure (2)
plot(ScopeData2.time,ScopeData2.signals(1,2).values,'LineWidth',1.5);hold on
ylabel('{\Delta}{\delta}')
hold on
PSS= 1; D=20;% Simulation with PSS
simopt = simset('solver','ode45','SrcWorkspace','Current','DstWorkspace','Current');  % Initialize sim options
[tout xout yout]=sim('SMIB_sim',[0 Time],simopt);  %#ok
plot(ScopeData2.time,ScopeData2.signals(1,2).values,'r','LineWidth',1.5);
legend('without PSS', 'with PSS')
ylabel('{\Delta}{\delta}')