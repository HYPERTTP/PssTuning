function [bestOrganism bestFitness]=RUNTHIS(ecosize,funnum)
clc
clear all
close all

% format short g
tic;
ecosize=20;
% Outputs: best organism/solution and best fitness
% Inputs: ecosystem/population size and # of benchmark problems
% Example: [A,B]=SOS (50,17), SOS will solve Sphere (F17) with 50 organisms
% (please see the "OBJECTIVE FUNCTIONS" and "SETUP" sub-functions)
%format compact
fprintf('-------------------------------------------------------------------------\n');
fprintf('  Symbiotic Organisms Search(SOS) for unconstrained benchmark problems\n');
fprintf('-------------------------------------------------------------------------\n\n');
funnum=1;
% --- Counters, Parameters & Matrix Initialization
[globalMin Lb1 Ub1 Lb2 Ub2 Lb3 Ub3 nd maxFE]=terminate(funnum);
fprintf(' Ecosystem Size: %d\t\tMaxFE: %d\t\tFunctionNumber: %d',ecosize,maxFE,funnum);
fprintf('\n\n');
fprintf('-------------------------------------------------------------------------\n\n');
FE=0;
% nd=9;                           % Function of Evaluation Counter
% eco=zeros(ecosize,nd);           % Ecosystem Matrix
fitness =zeros(ecosize,1);      % Fitness Matrix
% --- Ecosystem Initialization
eco=zeros(ecosize,nd);
iiv=1;
while (iiv~=ecosize+1)
    
    Kpss=rand*(Ub1-Lb1)+Lb1;
    T1=rand*(Ub2-Lb2)+Lb2;
    T2=rand*(Ub3-Lb3)+Lb3;
    eco(iiv,:)=horzcat(Kpss,T1,T2);
    
    %%% CALCULATION OF OBJECTIVE FUNCTION%%%%%%
    Re=0;
    Xe=0.4;
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
    D=20;
    ws=377;
    
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
    % =======================================
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
    
    %==================================
    
    x1=[-1/(K3*Tpdo)  -K4/Tpdo 0   1/Tpdo 0];
    x2=[ 0              0      ws    0    0];
    x3=[-K2/(2*H)  -K1/(2*H)   -D*ws/(2*H)  0 0];
    x4=[-KA*K6/TA  -KA*K5/TA  0 -1/TA  KA/TA];
    x5=[-K2*T1*Kpss/(2*H*T2)  -K1*T1*Kpss/(2*H*T2)  Kpss/T2  0 -1/T2];
    Apss = [ x1 ;x2 ; x3 ; x4; x5];
    eig0 = -1;
    eta0 = 0.2;
    [~,a2,a3 ] = damp(Apss);
    eigi= real (a3);
    
    dd1=length(a3);
    
    for no = 1:dd1
        j1 (no)= (eig0-eigi(no))^2;
        j2 (no)= (eta0-a2(no))^2;
        
    end
    
    F = 1*sum(j1)+ 1000*sum(j2);
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    fitness(iiv)=F;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  del_P_line
    %  del_FREQUENCY(4)
    %  %pauseclc
    iiv=iiv+1
end
% eco
% fitness
% pause

Pr=1;
while FE<maxFE
    for iiv=1:ecosize % Organisms' Looping
        
        % Update the best Organism
        [bestFitness,idx]=min(fitness);
        bestOrganism=eco(idx,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mutualism Phase
        % Choose organism j randomly other than organism i
        tv=1;
        while (tv~=2)
            jj=iiv;
            while iiv==jj
                seed=randperm(ecosize);
                jj=seed(1);
            end
            % Determine Mutual Vector & Beneficial Factor
            mutualVector=mean([eco(iiv,:);eco(jj,:)]);
            BF1=round(1+rand);
            % Calculate new solution after Mutualism Phase
            ecoNew1=eco(iiv,:)+rand(1,nd).*(bestOrganism-BF1.*mutualVector);
            clearvars Kpss T1 T2 Apss J1 J2 F
            Kpss=ecoNew1(1);
            T1= ecoNew1(2);
            T2= ecoNew1(3);
            
            if Kpss<=Lb1
                Kpss=Lb1;
            elseif Kpss>=Ub1
                Kpss=Ub1;
            end
            if T1<=Lb2
                T1=Lb2;
            elseif T1>=Ub2
                T1=Ub2;
            end
            if T2<=Lb3
                T2=Lb3;
            elseif T2>=Ub3
                T2=Ub3;
            end
            
            ecoNew1=horzcat(Kpss,T1,T2);
            x1=[-1/(K3*Tpdo)  -K4/Tpdo 0   1/Tpdo 0];
            x2=[ 0              0      ws    0    0];
            x3=[-K2/(2*H)  -K1/(2*H)   -D*ws/(2*H)  0 0];
            x4=[-KA*K6/TA  -KA*K5/TA  0 -1/TA  KA/TA];
            x5=[-K2*T1*Kpss/(2*H*T2)  -K1*T1*Kpss/(2*H*T2)  Kpss/T2  0 -1/T2];
            Apss = [ x1 ;x2 ; x3 ; x4; x5];
            eig0 = -1;
            eta0 = 0.2;
            [~,a2,a3 ] = damp(Apss);
            eigi= real (a3);
            dd1=length(a3);
            for no = 1:dd1
                j1 (no)= (eig0-eigi(no))^2;
                j2 (no)= (eta0-a2(no))^2;
            end
            
            F = 1*sum(j1)+ 1000*sum(j2);
            fitnessNew1(tv)=F;
            tv=tv+1;
        end
        %         fitnessNew1=realeig1;
%         ecoNew1
%         fitnessNew1
%         pause
        %%
        % for iiv=1:ecosize % Organisms' Looping
        
        % Update the best Organism
        %         [bestFitness,idx]=min(fitness); bestOrganism=eco(idx,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mutualism Phase
        % Choose organism j randomly other than organism i
        
        itt=1;
        while (itt~=2)
            jj=iiv;
            while iiv==jj
                seed=randperm(ecosize);
                jj=seed(1);
            end
            mutualVector=mean([eco(iiv,:);eco(jj,:)]);
            %             jj
            %             pause
            BF2=round(1+rand);
            ecoNew2=eco(jj,:)+rand(1,nd).*(bestOrganism-BF2.*mutualVector);
            clearvars Kpss T1 T2 Apss J1 J2 F
            Kpss=ecoNew2(1);
            T1= ecoNew2(2);
            T2= ecoNew2(3);
            
            if Kpss<=Lb1
                Kpss=Lb1;
            elseif Kpss>=Ub1
                Kpss=Ub1;
            end
            if T1<=Lb2
                T1=Lb2;
            elseif T1>=Ub2
                T1=Ub2;
            end
            if T2<=Lb3
                T2=Lb3;
            elseif T2>=Ub3
                T2=Ub3;
            end
            
            ecoNew2=horzcat(Kpss,T1,T2);
            
            x1=[-1/(K3*Tpdo)  -K4/Tpdo 0   1/Tpdo 0];
            x2=[ 0              0      ws    0    0];
            x3=[-K2/(2*H)  -K1/(2*H)   -D*ws/(2*H)  0 0];
            x4=[-KA*K6/TA  -KA*K5/TA  0 -1/TA  KA/TA];
            x5=[-K2*T1*Kpss/(2*H*T2)  -K1*T1*Kpss/(2*H*T2)  Kpss/T2  0 -1/T2];
            Apss = [ x1 ;x2 ; x3 ; x4; x5];
            eig0 = -1;
            eta0 = 0.2;
            Time = 20;
            dVref = 0.0;
            dTm= 0.20;
            tsim= 1.0005;
            [~,a2,a3 ] = damp(Apss);
            eigi= real (a3);
            
            dd1=length(a3);
            
            for no = 1:dd1
                j1 (no)= (eig0-eigi(no))^2;
                j2 (no)= (eta0-a2(no))^2;
                
            end
            
            F = 1*sum(j1)+ 1000*sum(j2);
            
            fitnessNew2(itt)=F;
            itt=itt+1;
        end
        %         fitnessNew2= realeig2;
%         ecoNew2
%         fitnessNew2
%         pause
        if fitnessNew1<fitness(iiv)
            fitness(iiv)=fitnessNew1;
            eco(iiv,:)=ecoNew1;
        end
        if fitnessNew2<fitness(jj)
            fitness(jj)=fitnessNew2;
            eco(jj,:)=ecoNew2;
        end
        % Increase the number of function evaluation counter
        FE=FE+2;
        % End of Mutualism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Commensialism Phase
        
        % Choose organism j randomly other than organism i
        it=1;
        while (it~=2)
            jj=iiv;
            while iiv==jj
                seed=randperm(ecosize);
                jj=seed(1);
            end
            ecoNew3=eco(iiv,:)+(rand(1,nd)*2-1).*(bestOrganism-eco(jj,:));
            clearvars Kpss T1 T2 Apss J1 J2 F
            Kpss=ecoNew3(1);
            T1= ecoNew3(2);
            T2= ecoNew3(3);
            
            if Kpss<=Lb1
                Kpss=Lb1;
            elseif Kpss>=Ub1
                Kpss=Ub1;
            end
            if T1<=Lb2
                T1=Lb2;
            elseif T1>=Ub2
                T1=Ub2;
            end
            if T2<=Lb3
                T2=Lb3;
            elseif T2>=Ub3
                T2=Ub3;
            end
            
            ecoNew2=horzcat(Kpss,T1,T2);
            
            x1=[-1/(K3*Tpdo)  -K4/Tpdo 0   1/Tpdo 0];
            x2=[ 0              0      ws    0    0];
            x3=[-K2/(2*H)  -K1/(2*H)   -D*ws/(2*H)  0 0];
            x4=[-KA*K6/TA  -KA*K5/TA  0 -1/TA  KA/TA];
            x5=[-K2*T1*Kpss/(2*H*T2)  -K1*T1*Kpss/(2*H*T2)  Kpss/T2  0 -1/T2];
            Apss = [ x1 ;x2 ; x3 ; x4; x5];
            eig0 = -1;
            eta0 = 0.2;
            Time = 20;
            dVref = 0.0;
            dTm= 0.20;
            tsim= 1.0005;
            [~,a2,a3 ] = damp(Apss);
            eigi= real (a3);
            
            dd1=length(a3);
            
            for no = 1:dd1
                j1 (no)= (eig0-eigi(no))^2;
                j2 (no)= (eta0-a2(no))^2;
                
            end
            
            F = 1*sum(j1)+ 1000*sum(j2);
            fitnessNew3(it)=F;
            it=it+1;
        end
        
%         ecoNew3
%         fitnessNew3
%         pause
        %
        % Evaluate the fitness of the new solution
        
        
        % Accept the new solution if the fitness is better
        if fitnessNew3<fitness(iiv)
            fitness(iiv)=fitnessNew3;
            eco(iiv,:)=ecoNew3;
        end
        % Increase the number of function evaluation counter
        FE=FE+1;
        %     end
        % End of Commensalism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parasitism Phase
        
        % Choose organism j randomly other than organism i
        ittt=1;
        while ittt~=2
            jj=iiv;
            while iiv==jj
                seed=randperm(ecosize);
                jj=seed(1);
            end
            % Determine Parasite Vector & Calculate the fitness
            parasiteVector=eco(iiv,:);
            seed=randperm(nd);
            pick=seed(1:ceil(rand*nd));  % select random dimension
            ub=[Ub1 Ub2 Ub3];
            lb=[Lb1 Lb2 Lb3];
            parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
            clearvars Kpss T1 T2 Apss J1 J2 F
            Kpss=parasiteVector(1);
            T1= parasiteVector(2);
            T2= parasiteVector(3);
            
            if Kpss<=Lb1
                Kpss=Lb1;
            elseif Kpss>=Ub1
                Kpss=Ub1;
            end
            if T1<=Lb2
                T1=Lb2;
            elseif T1>=Ub2
                T1=Ub2;
            end
            if T2<=Lb3
                T2=Lb3;
            elseif T2>=Ub3
                T2=Ub3;
            end
            
            parasiteVector=horzcat(Kpss,T1,T2);
            
            x1=[-1/(K3*Tpdo)  -K4/Tpdo 0   1/Tpdo 0];
            x2=[ 0              0      ws    0    0];
            x3=[-K2/(2*H)  -K1/(2*H)   -D*ws/(2*H)  0 0];
            x4=[-KA*K6/TA  -KA*K5/TA  0 -1/TA  KA/TA];
            x5=[-K2*T1*Kpss/(2*H*T2)  -K1*T1*Kpss/(2*H*T2)  Kpss/T2  0 -1/T2];
            Apss = [ x1 ;x2 ; x3 ; x4; x5];
            eig0 = -1;
            eta0 = 0.2;
            Time = 20;
            dVref = 0.0;
            dTm= 0.20;
            tsim= 1.0005;
            [~,a2,a3 ] = damp(Apss);
            eigi= real (a3);
            
            dd1=length(a3);
            
            for no = 1:dd1
                j1 (no)= (eig0-eigi(no))^2;
                j2 (no)= (eta0-a2(no))^2;
                
            end
            
            F = 1*sum(j1)+ 1000*sum(j2);
            fitnessParasite(ittt)=F;
            ittt=ittt+1;
        end
%         parasiteVector
%         fitnessParasite
%         pause
        if fitnessParasite < fitness(jj)
            fitness(jj)=fitnessParasite;
            eco(jj,:)=parasiteVector;
        end
        % Increase the number of function evaluation counter
        FE=FE+1
    end
    % End of Parasitism Phase
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     end % End of Organisms' Looping
    % Checking the termination criteria
    %     if bestFitness<globalMin
    %         break
    %     end
    FE
    bestFitness
    bestOrganism
    %     pause
    ggp(Pr,:)=bestFitness;
    hhp1(Pr,:)=Pr;
    Pr=Pr+1;
end
% End of Main Looping
% --- Update the best Organism
[bestFitness,idx]=min(fitness)
bestOrganism=eco(idx,:)
figure(3)
plot(hhp1,ggp)
% --- Display the result
%disp(['Funnum: ', num2str(funnum)])
disp(['Best Fitness: ', num2str(bestFitness)])
disp(['Best Organism ', num2str(bestOrganism)])
toc;
% pause
end
% function fitness=fobj(del_FREQUENCY,DELPTIE_new) % GO TO LINE NO. 587
% ecosize=5;
% for ii=1:ecosize
%    err(ii)=del_FREQUENCY(4)^2+del_FREQUENCY(6)^2+del_FREQUENCY(10)^2+DELPTIE_new(1)^2+DELPTIE_new(2)^2+DELPTIE_new(3)^2;
% end
% fitness=err;
% end

function [globalMin Lb1 Ub1 Lb2 Ub2 Lb3 Ub3 nd maxFE]=terminate(funnum)
maxFE=4000; % maximum number of function evaluation
Tol_opti =1e-12;
Lb1=[0.1];       % TCSC
Ub1=[50];
Lb2=[0.2];       % Ts
Ub2=[1.5];
Lb3=[0.02];       % Tt
Ub3=[0.15];

nd=3;                       % No. of parameters to be optimized ( )
globalMin = 0;
globalMin = globalMin+Tol_opti;
end

