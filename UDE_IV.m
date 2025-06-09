function [Sol1, Sol2, Sol3,Save_data] = UDE_IV(parameters)

% Extracted Parameters %%%%%%
run = parameters.run;
D = parameters.D;
func_No = parameters.func_No;
Np = parameters.Np;
x_Max = parameters.x_Max;
x_Min = parameters.x_Min;
frac_T = parameters.frac_T;
frac_T1 = parameters.frac_T1;
Lp = parameters.Lp;
max_FEs = parameters.max_FEs;
H = parameters.H;
memory = parameters.memory;
plot_flag = parameters.plot_flag;
para_adap = parameters.para_adap;
sort_flag = parameters.sort_flag;
Max_gen = parameters.Max_gen;
pc = parameters.pc;
numcon = parameters.numcon;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization %%%%%%%%%%%%
Pop = x_Min*ones(Np,D) + (x_Max-x_Min)*rand(Np,D);
[Func_Val, Cons_Viol] = get_Fitness(Pop, func_No);
FEs = 150;

% Sorting the initial population according to SOF, leads to migration of solutions between sub-populations %
sorted_Id = sort1(Func_Val, Cons_Viol);
Pop = Pop(sorted_Id, :);
Func_Val = Func_Val(sorted_Id);
Cons_Viol = Cons_Viol(sorted_Id);

% Parameter Pool of CoDE %%%%%%%%%%%%
%F_Pool =  [  1,   1, 0.8];
%CR_Pool = [0.1, 0.9, 0.2];

% Parameter Pool of UDE 2017 Conference %%%%%%%%%%%%
%F_Pool =  [0.9, 0.5];
%CR_Pool = [0.9, 0.5];
feasibility_Percentage = [];
Constraint_violation=[];
if para_adap == 1
    M_CR1 = 0.5*ones(1,H); M_CR3 = 0.5*ones(1,H);
    CR = ones(3,Np);
    M_F1 = 0.5*ones(1,H); M_F2 = 0.5*ones(1,H); M_F3 = 0.5*ones(1,H);
    F = ones(3,Np);
end

k1 = 1; k2 = 1; k3 = 1;

Gen = 1;
SProp = 0.40;
SG = 25;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Temp_Pop = zeros(Np,D);
Temp_Val = zeros(Np,1);
Temp_Cons = zeros(Np,1);
interval_FEs = [2000*D, 10000*D, 20000*D];
interval_Count = 1;

best_Vec = Pop(1,:);
best_Func_Val = Func_Val(1);
best_Cons_Viol = Cons_Viol(1);

Save_data.Best_Func_Val=[];
Save_data.Best_Cons_Viol=[];

NW = zeros(Max_gen+1,3);
NL = zeros(Max_gen+1,3);
Temp_Id = zeros(Max_gen+1,3);
CVmax = max(Cons_Viol);
epsilon0 = CVmax;
cp=(-log(epsilon0)-6)/log(1-pc);
if(cp > 33)
    cp = 33;
end
Archive=[];
index=zeros(Max_gen,Np);

Count_track_No_of_Stagnated_Vectors = zeros(Max_gen,1);

XG = zeros(Max_gen,D);
dG = zeros(Max_gen,1);

XG(Gen,:) = (1/Np)*sum(Pop);
dist = 0;
for i = 1 : Np
    dist = dist + norm(Pop(i,:)-XG(Gen,:));        
end
dG(Gen,1) = (1/Np)*dist;

Archive_TV=[];Archive_Val_TV=[];Archive_Cons_TV=[];
while FEs <= max_FEs
 %% Initialize temporary archive in each generation
  Temp_A=[];Temp_A_Val=[];Temp_A_Cons=[]; 
  
    Gen = Gen + 1;
    if( Gen/Max_gen <= pc)
        epsilon = epsilon0*((1 - (Gen/Max_gen))^cp);
        %epsilonP(Gen) = epsilon;
    else
        epsilon = 0;
        %epsilonP(Gen) = epsilon;
    end
    
    prob = (Np-1:-1:0)/Np;
    
    % For Top T Population %%%%%%%%%%%
    T = frac_T*Np;
    T1 = frac_T1*Np;
    S1_Pop = zeros(T, D);
    S2_Pop = zeros(T, D);
    S3_Pop = zeros(T, D);
    
    for i = 1:T1
        
        % DE/rand/1 %%%%%%%%%%%%%%%%%%%
        r1 = randi(Np);
        while r1 == i
            r1 = randi(Np);
        end
        
        r2 = randi(Np);
        while r2 == r1 || r2 == i
            r2 = randi(Np);
        end
        
        r3 = randi(Np);
        while r3 == r1 || r3 == r2 || r3 == i
            r3 = randi(Np);
        end
        
        if para_adap == 1
            temp_rand = randi(H);
            CR(1,i) = normrnd(M_CR1(temp_rand), 0.1);
            CR(1,i) = min(CR(1,i), 1); CR(1,i) = max(CR(1,i), 0);
            
            F(1,i) = cauchyrand(M_F1(temp_rand), 0.1);
            while F(1,i)<=0
                F(1,i) = cauchyrand(M_F1(temp_rand), 0.1);
            end
            F(1,i) = min(F(1,i), 1);
        else
            temp_rand = randi(length(F_Pool));
            F(1,i) = F_Pool(temp_rand);
            CR(1,i) = CR_Pool(temp_rand);
        end
        
        % Mutation
        S1_Pop(i,:) = Pop(r1,:) + F(1,i)*(Pop(r2,:)-Pop(r3,:));
        
        % Crossover
        j_rand = randi(D);
        temp = rand(1, D);
        U1 = zeros(1,D);
        U1(temp <= CR(1,i)) = S1_Pop(i,temp <= CR(1,i));
        U1(j_rand) = S1_Pop(i,j_rand);
        U1(temp > CR(1,i)) = Pop(i, temp > CR(1,i));
        
        S1_Pop(i,:) = U1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % DE/current-to-rand/1 %%%%%%%%
        r1 = randi(Np);
        while r1 == i
            r1 = randi(Np);
        end
        
        r2 = randi(Np);
        while r2 == r1 || r2 == i
            r2 = randi(Np);
        end
        
        r3 = randi(Np);
        while r3 == r1 || r3 == r2 || r3 == i
            r3 = randi(Np);
        end
        
        if para_adap == 1
            temp_rand = randi(H);
            F(2,i) = cauchyrand(M_F2(temp_rand), 0.1);
            while F(2,i)<=0
                F(2,i) = cauchyrand(M_F2(temp_rand), 0.1);
            end
            F(2,i) = min(F(2,i), 1);
        else
            temp_rand = randi(length(F_Pool));
            F(2,i) = F_Pool(temp_rand);
            CR(2,i) = CR_Pool(temp_rand);
        end
        
        % Mutation
        S2_Pop(i,:) = Pop(i,:) + F(2,i)*(Pop(r1,:)-Pop(i,:)) + F(2,i)*(Pop(r2,:)-Pop(r3,:));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % DE/current-to-pbest/1 %%%%%%%
        a = 5; b = 20;
        p = round((b-a)*rand + a);
        %p = 5;
        
        r1 = randi(p);
        while r1 == i
            r1 = randi(p);
        end
        
        % r1 = randi(Np);
        % while rand >= prob(r1) || r1 == i
        % 	r1 = randi(Np);
        % end
        
        r2 = randi(Np);
        while r2 == r1 || r2 == i
            r2 = randi(Np);
        end
        
        r3 = randi(Np);
        while r3 == r1 || r3 == r2 || r3 == i
            r3 = randi(Np);
        end
        
        if para_adap == 1
            temp_rand = randi(H);
            CR(3,i) = normrnd(M_CR3(temp_rand), 0.1);
            CR(3,i) = min(CR(3,i), 1); CR(3,i) = max(CR(3,i), 0);
            
            F(3,i) = cauchyrand(M_F3(temp_rand), 0.1);
            while F(3,i)<=0
                F(3,i) = cauchyrand(M_F3(temp_rand), 0.1);
            end
            F(3,i) = min(F(3,i), 1);
        else
            temp_rand = randi(length(F_Pool));
            F(3,i) = F_Pool(temp_rand);
            CR(3,i) = CR_Pool(temp_rand);
        end
        
        % Mutation
        S3_Pop(i,:) = Pop(i,:) + F(3,i)*(Pop(r1,:)-Pop(i,:)) + F(3,i)*(Pop(r2,:)-Pop(r3,:));
        
        % Crossover
        j_rand = randi(D);
        temp = rand(1, D);
        U1 = zeros(1,D);
        U1(temp <= CR(3,i)) = S3_Pop(i,temp <= CR(3,i));
        U1(j_rand) = S3_Pop(i,j_rand);
        U1(temp > CR(3,i)) = Pop(i, temp > CR(3,i));
        
        S3_Pop(i,:) = U1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    S1_Pop(S1_Pop < x_Min) = x_Min;
    S1_Pop(S1_Pop > x_Max) = x_Max;
    
    S2_Pop(S2_Pop < x_Min) = x_Min;
    S2_Pop(S2_Pop > x_Max) = x_Max;
    
    S3_Pop(S3_Pop < x_Min) = x_Min;
    S3_Pop(S3_Pop > x_Max) = x_Max;
    
    Trial_Vec = [S1_Pop; S2_Pop; S3_Pop];
    [Trial_Val, Trial_Cons] = get_Fitness(Trial_Vec, func_No);
    FEs = FEs + 3*T1;
    
    if para_adap == 1
        SCR1 = zeros(T1,1); SCR3 = zeros(T1,1);
        SF1 = zeros(T1,1); SF2 = zeros(T1,1); SF3 = zeros(T1,1);
        W1 = zeros(T1,1); W2 = zeros(T1,1); W3 = zeros(T1,1);
        id1 = 0; id2 = 0; id3 = 0;
    end
    
    for i = 1:T1
        temp_Id = [i, T1+i, 2*T1+i];
        
        Val = Trial_Val(temp_Id);
        Cons = Trial_Cons(temp_Id);
        
        Sorted_Id = sort1(Val, Cons);
        min_Id = Sorted_Id(1);
        
        NW(Gen, min_Id) = NW(Gen, min_Id)+1;
        
        Temp_Pop(i,:) = Trial_Vec(temp_Id(min_Id),:);
        
        Temp_Val(i) = Trial_Val(temp_Id(min_Id));
        Temp_Cons(i) = Trial_Cons(temp_Id(min_Id));
        Temp_Id(Gen, i) = min_Id;
        
        if para_adap == 1
            
            if(Trial_Cons(temp_Id(1)) > 0 && Cons_Viol(i) > 0)
                if(Trial_Cons(temp_Id(1)) < Cons_Viol(i))
                    id1 = id1 + 1;
                    SCR1(id1) = CR(1,i);
                    SF1(id1) = F(1,i);
                    W1(id1) = (Cons_Viol(i)/numcon) - (Trial_Cons(temp_Id(1))/numcon);
                end
            end
            
            if(Trial_Cons(temp_Id(1)) == 0 && Cons_Viol(i) > 0)
                id1 = id1 + 1;
                SCR1(id1) = CR(1,i);
                SF1(id1) = F(1,i);
                W1(id1) = (Cons_Viol(i)/numcon) - 0;
            end
            
            if(Trial_Cons(temp_Id(1)) == 0 && Cons_Viol(i) == 0)
                if(Trial_Val(temp_Id(1)) < Func_Val(i))
                    id1 = id1 + 1;
                    SCR1(id1) = CR(1,i);
                    SF1(id1) = F(1,i);
                    W1(id1) = Func_Val(i)-Trial_Val(temp_Id(1));
                end
            end
            
            
            if(Trial_Cons(temp_Id(2)) > 0 && Cons_Viol(i) > 0)
                if(Trial_Cons(temp_Id(2)) < Cons_Viol(i))
                    id2 = id2 + 1;
                    SF2(id2) = F(2,i);
                    W2(id2) = (Cons_Viol(i)/numcon) - (Trial_Cons(temp_Id(2))/numcon);
                end
            end
            
            if(Trial_Cons(temp_Id(2)) == 0 && Cons_Viol(i) > 0)
                id2 = id2 + 1;
                SF2(id2) = F(2,i);
                W2(id2) = (Cons_Viol(i)/numcon) - 0;
            end
            
            if(Trial_Cons(temp_Id(2)) == 0 && Cons_Viol(i) == 0)
                if(Trial_Val(temp_Id(2)) < Func_Val(i))
                    id2 = id2 + 1;
                    SF2(id2) = F(2,i);
                    W2(id2) = Func_Val(i)-Trial_Val(temp_Id(2));
                end
            end
            
            if(Trial_Cons(temp_Id(3)) > 0 && Cons_Viol(i) > 0)
                if(Trial_Cons(temp_Id(3)) < Cons_Viol(i))
                    id3 = id3 + 1;
                    SCR3(id3) = CR(3,i);
                    SF3(id3) = F(3,i);
                    W3(id3) = (Cons_Viol(i)/numcon) - (Trial_Cons(temp_Id(3))/numcon);
                end
            end
            
            if(Trial_Cons(temp_Id(3)) == 0 && Cons_Viol(i) > 0)
                id3 = id3 + 1;
                SCR3(id3) = CR(3,i);
                SF3(id3) = F(3,i);
                W3(id3) = (Cons_Viol(i)/numcon) - 0;
            end
            
            if(Trial_Cons(temp_Id(3)) == 0 && Cons_Viol(i) == 0)
                if(Trial_Val(temp_Id(3)) < Func_Val(i))
                    id3 = id3 + 1;
                    SCR3(id3) = CR(3,i);
                    SF3(id3) = F(3,i);
                    W3(id3) = Func_Val(i)-Trial_Val(temp_Id(3));
                end
            end
            
        end
    end
    
    
    if para_adap == 1
        SCR1 = SCR1(1:id1); SF1 = SF1(1:id1); W1 = W1(1:id1);
        SF2 = SF2(1:id2); W2 = W2(1:id2);
        SCR3 = SCR3(1:id3); SF3 = SF3(1:id3); W3 = W3(1:id3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % For Remaining Polulation %%%%%%%%
    if Gen < (Lp + 1)
        Gen;
        %SR = 0.33*ones(1,3);
        SR(Gen,1) = 0.33;
        SR(Gen,2) = 0.33;
        SR(Gen,3) = 0.33;
    elseif(Gen == (Lp+1))
        NW1 = 0; NW2 = 0; NW3 = 0;
        NL1 = 0; NL2 = 0; NL3 = 0;
        %NR1 = 0; NR2 = 0; NR3 = 0;
        for g = (Gen - Lp) : 1: (Gen - 1)
            NW1 = NW1 + NW(g,1);
            NW2 = NW2 + NW(g,2);
            NW3 = NW3 + NW(g,3);
            NL1 = NL1 + NL(g,1);
            NL2 = NL2 + NL(g,2);
            NL3 = NL3 + NL(g,3);
        end
        NR1 = (NW1/(NW1+NL1)) + 0.01;
        NR2 = (NW2/(NW2+NL2)) + 0.01;
        NR3 = (NW3/(NW3+NL3)) + 0.01;
        SR(Gen,1) = NR1/(NR1 + NR2 + NR3);
        SR(Gen,2) = NR2/(NR1 + NR2 + NR3);
        SR(Gen,3) = NR3/(NR1 + NR2 + NR3);
    elseif(Gen > (Lp+1))
        NW1 = NW1 - NW(Gen - Lp-1,1) + NW(Gen-1,1);
        NW2 = NW2 - NW(Gen - Lp-1,2) + NW(Gen-1,2);
        NW3 = NW3 - NW(Gen - Lp-1,3) + NW(Gen-1,3);
        NL1 = NL1 - NL(Gen - Lp-1,1) + NL(Gen-1,1);
        NL2 = NL2 - NL(Gen - Lp-1,2) + NL(Gen-1,2);
        NL3 = NL3 - NL(Gen - Lp-1,3) + NL(Gen-1,3);
        NR1 = (NW1/(NW1+NL1)) + 0.01;
        NR2 = (NW2/(NW2+NL2)) + 0.01;
        NR3 = (NW3/(NW3+NL3)) + 0.01;
        SR(Gen,1) = NR1/(NR1 + NR2 + NR3);
        SR(Gen,2) = NR2/(NR1 + NR2 + NR3);
        SR(Gen,3) = NR3/(NR1 + NR2 + NR3);
    end
    
    for i = T1+1:Np
        temp_rand = rand;
        
        if temp_rand <= SR(Gen,1)
            Temp_Id(Gen,i) = 1;
            % DE/rand/1 %%%%%%%%%%%%%%%
            r1 = randi(Np);
            while r1 == i
                r1 = randi(Np);
            end
            
            r2 = randi(Np);
            while r2 == r1 || r2 == i
                r2 = randi(Np);
            end
            
            r3 = randi(Np);
            while r3 == r1 || r3 == r2 || r3 == i
                r3 = randi(Np);
            end
            
            if para_adap == 1
                temp_rand = randi(H);
                CR(1,i) = normrnd(M_CR1(temp_rand), 0.1);
                CR(1,i) = min(CR(1,i), 1); CR(1,i) = max(CR(1,i), 0);
                
                F(1,i) = cauchyrand(M_F1(temp_rand), 0.1);
                while F(1,i)<=0
                    F(1,i) = cauchyrand(M_F1(temp_rand), 0.1);
                end
                F(1,i) = min(F(1,i), 1);
            else
                temp_rand = randi(length(F_Pool));
                F(1,i) = F_Pool(temp_rand);
                CR(1,i) = CR_Pool(temp_rand);
            end
            
            % Mutation
            Trial_Vec = Pop(r1,:) + F(1,i)*(Pop(r2,:)-Pop(r3,:));
            
            % Crossover
            j_rand = randi(D);
            temp = rand(1, D);
            U1 = zeros(1,D);
            U1(temp <= CR(1,i)) = Trial_Vec(temp <= CR(1,i));
            U1(j_rand) = Trial_Vec(j_rand);
            U1(temp > CR(1,i)) = Pop(i, temp > CR(1,i));
            
            Trial_Vec = U1;
            
        elseif temp_rand <= (SR(Gen,1)+SR(Gen,2))
            Temp_Id(Gen,i) = 2;
            % DE/current-to-rand/1 %%%%
            r1 = randi(Np);
            while r1 == i
                r1 = randi(Np);
            end
            
            r2 = randi(Np);
            while r2 == r1 || r2 == i
                r2 = randi(Np);
            end
            
            r3 = randi(Np);
            while r3 == r1 || r3 == r2 || r3 == i
                r3 = randi(Np);
            end
            
            if para_adap == 1
                temp_rand = randi(H);
                F(2,i) = cauchyrand(M_F2(temp_rand), 0.1);
                while F(2,i)<=0
                    F(2,i) = cauchyrand(M_F2(temp_rand), 0.1);
                end
                F(2,i) = min(F(2,i), 1);
            else
                temp_rand = randi(length(F_Pool));
                F(2,i) = F_Pool(temp_rand);
                CR(2,i) = CR_Pool(temp_rand);
            end
            
            % Mutation
            Trial_Vec = Pop(i,:) + F(2,i)*(Pop(r1,:)-Pop(i,:)) + F(2,i)*(Pop(r2,:)-Pop(r3,:));
            
        else
            Temp_Id(Gen,i) = 3;
            % DE/current-to-pbest/1 %%%%%%%
            a = 5; b = 20;
            p = round((b-a)*rand + a);
            %p = 5;
            
            r1 = randi(p);
            
            while r1 == i
                r1 = randi(p);
            end
            
            r0 = randi(Np);
            while r0 == r1 || r0 == i
                r0 = randi(Np);
            end
            
            r2 = randi(Np);
            while r2==r0 || r2 == r1 || r2 == i
                r2 = randi(Np);
            end
            
            r3 = randi(Np);
            while r3==r0 || r3 == r1 || r3 == r2 || r3 == i
                r3 = randi(Np);
            end
            
            if para_adap == 1
                temp_rand = randi(H);
                CR(3,i) = normrnd(M_CR3(temp_rand), 0.1);
                CR(3,i) = min(CR(3,i), 1); CR(3,i) = max(CR(3,i), 0);
                
                F(3,i) = cauchyrand(M_F3(temp_rand), 0.1);
                while F(3,i)<=0
                    F(3,i) = cauchyrand(M_F3(temp_rand), 0.1);
                end
                F(3,i) = min(F(3,i), 1);
            else
                temp_rand = randi(length(F_Pool));
                F(3,i) = F_Pool(temp_rand);
                CR(3,i) = CR_Pool(temp_rand);
            end
            
          % Hybrid Mutations
            if rand<0.75
                % DE/current-to-pbest/1
                Trial_Vec = Pop(i,:) + F(3,i)*(Pop(r1,:)-Pop(i,:)) + F(3,i)*(Pop(r2,:)-Pop(r3,:));
            else
                % DE/rand-to-pbest/1
                Trial_Vec = Pop(r0,:) + F(3,i)*(Pop(r1,:)-Pop(r0,:)) + F(3,i)*(Pop(r2,:)-Pop(r3,:));
            end
            
            % Crossover
            j_rand = randi(D);
            temp = rand(1, D);
            U1 = zeros(1,D);
            U1(temp <= CR(3,i)) = Trial_Vec(temp <= CR(3,i));
            U1(j_rand) = Trial_Vec(j_rand);
            U1(temp > CR(3,i)) = Pop(i, temp > CR(3,i));
            
            Trial_Vec = U1;
        end
        
        Trial_Vec(Trial_Vec < x_Min) = x_Min;
        Trial_Vec(Trial_Vec > x_Max) = x_Max;
        
        Temp_Pop(i,:) = Trial_Vec;
        
    end
    
    [Temp_Val(T1+1:Np), Temp_Cons(T1+1:Np)] = get_Fitness(Temp_Pop(T1+1:Np, :), func_No);
    FEs = FEs + (Np-T1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % One-to-One Comparison (Replacement/Selection)  %%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:Np
        %%% Replacement if strategy was DE/rand/1 and recording its win and
        %%% losses in creating successful offspring
        if(Temp_Id(Gen,i)==1)
            if( epsilon > 0 && Temp_Cons(i) <= epsilon && Cons_Viol(i) <= epsilon )
                if( Temp_Val(i) < Func_Val(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,1) = NW(Gen,1) + 1;
                    index(Gen,i)=0;
                else
                    NL(Gen,1) = NL(Gen,1) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            elseif( Temp_Cons(i) == Cons_Viol(i) )
                if( Temp_Val(i) < Func_Val(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,1) = NW(Gen,1) + 1;
                    index(Gen,i)=0;
                else
                    NL(Gen,1) = NL(Gen,1) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            else
                if( Temp_Cons(i) < Cons_Viol(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,1) = NW(Gen,1) + 1;
                    index(Gen,i)=0;
                elseif( Temp_Cons(i) > Cons_Viol(i) )
                    NL(Gen,1) = NL(Gen,1) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            end
        end
        
        %%% Replacement if strategy was DE/current-to-rand/1 and recording its win and
        %%% losses in creating successful offspring
        if(Temp_Id(Gen,i)==2)
            if( epsilon > 0 && Temp_Cons(i) <= epsilon && Cons_Viol(i) <= epsilon )
                if( Temp_Val(i) < Func_Val(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,2) = NW(Gen,2) + 1;
                    index(Gen,i)=0;
                else
                    NL(Gen,2) = NL(Gen,2) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            elseif( Temp_Cons(i) == Cons_Viol(i) )
                if( Temp_Val(i) < Func_Val(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,2) = NW(Gen,2) + 1;
                    index(Gen,i)=0;
                else
                    NL(Gen,2) = NL(Gen,2) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            else
                if( Temp_Cons(i) < Cons_Viol(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,2) = NW(Gen,2) + 1;
                    index(Gen,i)=0;
                elseif( Temp_Cons(i) > Cons_Viol(i) )
                    NL(Gen,2) = NL(Gen,2) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            end
        end
        
        %%% Replacement if strategy was DE/current-to-pbest/1 and recording its win and
        %%% losses in creating successful offspring
        if(Temp_Id(Gen,i)==3)
            if( epsilon > 0 && Temp_Cons(i) <= epsilon && Cons_Viol(i) <= epsilon )
                if( Temp_Val(i) < Func_Val(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,3) = NW(Gen,3) + 1;
                    index(Gen,i)=0;
                else
                    NL(Gen,3) = NL(Gen,3) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            elseif( Temp_Cons(i) == Cons_Viol(i) )
                if( Temp_Val(i) < Func_Val(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,3) = NW(Gen,3) + 1;
                    index(Gen,i)=0;
                else
                    NL(Gen,3) = NL(Gen,3) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            else
                if( Temp_Cons(i) < Cons_Viol(i) )
                    if ((Gen >= (SG+1)) && (index(Gen-1,i)>SG))
                        Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
                    end
                    Pop(i,:) = Temp_Pop(i,:);
                    Cons_Viol(i) = Temp_Cons(i);
                    Func_Val(i) = Temp_Val(i);
                    NW(Gen,3) = NW(Gen,3) + 1;
                    index(Gen,i)=0;
                elseif( Temp_Cons(i) > Cons_Viol(i) )
                    NL(Gen,3) = NL(Gen,3) + 1;
                    index(Gen,i)=index(Gen-1,i)+1;
                    Temp_A=[Temp_A;Temp_Pop(i,:)];
                    Temp_A_Val=[Temp_A_Val;Temp_Val(i)];
                    Temp_A_Cons=[Temp_A_Cons;Temp_Cons(i)];
                end
            end
        end
        
    end
    
    sorted_Id = sort1(Func_Val, Cons_Viol);
    Pop = Pop(sorted_Id, :);
    Func_Val = Func_Val(sorted_Id);
    Cons_Viol = Cons_Viol(sorted_Id);
    index(Gen,:) = index(Gen,sorted_Id);
    tt(Gen)=size(Temp_A,1);
    
    %% Store temporary Archive values to the Main Archive
    Archive_TV=[Archive_TV;Temp_A];
    Archive_Val_TV =[Archive_Val_TV;Temp_A_Val];
    Archive_Cons_TV =[Archive_Cons_TV;Temp_A_Cons];

    %% Turncated Archive size (Archive size will always =Np)
    if size(Archive_TV,1)>Np
        Sorted_Id_A=sort1(Archive_Val_TV,Archive_Cons_TV);
        Archive_TV =Archive_TV(Sorted_Id_A,:);
        
        Archive_Val_TV =Archive_Val_TV(Sorted_Id_A);
        Archive_Cons_TV =Archive_Cons_TV(Sorted_Id_A);
        Archive_TV(Np+1:end,:)=[];
        Archive_Val_TV(Np+1:end)=[];
        Archive_Cons_TV(Np+1:end)=[];
    end
    Archive_TV;

    % Finding centroid of the population
    if(Gen <= Max_gen)
        XG(Gen,:) = (1/Np)*sum(Pop);
        dist = 0;
        for i = 1 : Np
            dist = dist + norm(Pop(i,:)-XG(Gen,:));        
        end
        dG(Gen,1) = (1/Np)*dist;
    end

    if(Gen <= Max_gen)
        for i=1:Np
            if index(Gen,i)>SG
                Count_track_No_of_Stagnated_Vectors(Gen,1) = Count_track_No_of_Stagnated_Vectors(Gen,1)+1;
            end
        end
    end


if ((Gen <= Max_gen) && (Count_track_No_of_Stagnated_Vectors(Gen,1) > (SProp)*Np) && (dG(Gen,1)>1))
    for i=1:Np
        if index(Gen,i)> SG %SG=Lp
            id=unidrnd(size(Archive_TV,1));
            Pop(i,:)= Archive_TV(id,:);
            Func_Val(i)=Archive_Val_TV(id);
            Cons_Viol(i)=Archive_Cons_TV(id);
            index(Gen,i) = 0;
            break
        end
    end
end

    % Finding the best according to SOF in the current population
    cur_best_Vec_Id = sort1(Func_Val, Cons_Viol);
    cur_best_Vec_Id = cur_best_Vec_Id(1);
    cur_Best_Vec = Pop(cur_best_Vec_Id, :);
    cur_Best_Func_Val = Func_Val(cur_best_Vec_Id);
    cur_Best_Cons_Viol = Cons_Viol(cur_best_Vec_Id);
    
    % Archiving the best solution found so far based on SOF %
    if memory == 1
        if (cur_Best_Cons_Viol < best_Cons_Viol) || (cur_Best_Cons_Viol == best_Cons_Viol && cur_Best_Func_Val < best_Func_Val)
            best_Vec = cur_Best_Vec;
            best_Func_Val = cur_Best_Func_Val;
            best_Cons_Viol = cur_Best_Cons_Viol;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Updation of MCR and MF %%%%%%%%%%
    if para_adap == 1
        %if isempty(SCR1)==0 && isempty(SF1)==0
        if id1 > 0
            k1;
            M_CR1(k1) = update_MCR(SCR1, W1);
            M_F1(k1) = update_MF(SF1, W1);
            k1 = mod(k1,H) + 1;
        end
        
        %if isempty(SF2)==0
        if id2 > 0
            k2;
            M_F2(k2) = update_MF(SF2, W2);
            k2 = mod(k2,H) + 1;
        end
        
        %if isempty(SF3)==0
        if id3 > 0
            k3;
            M_CR3(k3) = update_MCR(SCR3, W3);
            M_F3(k3) = update_MF(SF3, W3);
            k3 = mod(k3,H) + 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Storing the interval results
    if FEs >= interval_FEs(1) && interval_Count == 1
        if memory == 1
            Sol1.vector = best_Vec;
            Sol1.func_Val = best_Func_Val;
            Sol1.cons_Viol = best_Cons_Viol;
        else
            Sol1.vector = cur_Best_Vec;
            Sol1.func_Val = cur_Best_Func_Val;
            Sol1.cons_Viol = cur_Best_Cons_Viol;
        end
        interval_Count = interval_Count + 1;
    elseif FEs >= interval_FEs(2) && interval_Count == 2
        if memory == 1
            Sol2.vector = best_Vec;
            Sol2.func_Val = best_Func_Val;
            Sol2.cons_Viol = best_Cons_Viol;
        else
            Sol2.vector = cur_Best_Vec;
            Sol2.func_Val = cur_Best_Func_Val;
            Sol2.cons_Viol = cur_Best_Cons_Viol;
        end
        interval_Count = interval_Count + 1;
    elseif FEs >= interval_FEs(3)
        if memory == 1
            Sol3.vector = best_Vec;
            Sol3.func_Val = best_Func_Val;
            Sol3.cons_Viol = best_Cons_Viol;
        else
            Sol3.vector = cur_Best_Vec;
            Sol3.func_Val = cur_Best_Func_Val;
            Sol3.cons_Viol = cur_Best_Cons_Viol;
        end
    end
    
    if(Gen == Max_gen)
        for i=1:Np
            if index(Gen,i)>SG
                Archive=[Archive;i, Gen - index(Gen-1,i),Gen, index(Gen-1,i), Func_Val(i), Cons_Viol(i)];
            end
        end
    end
    
    	if Gen==1 || mod(Gen,25)==0
    		if memory == 1
    			fprintf('%2d.%2d| Gen: %d, Func Val: %f, CV: %f\n', func_No, run, Gen, best_Func_Val, best_Cons_Viol);
    			elseif memory == 0
    			fprintf('%2d.%2d| Gen: %d, Func Val: %f, CV: %f\n', func_No, run, Gen, cur_Best_Func_Val, cur_Best_Cons_Viol);
    		end
        end
    
    Best_Cons_Viol(Gen) = best_Cons_Viol;
    Best_Func_Val(Gen) = best_Func_Val;
    feasibility_Percentage = [feasibility_Percentage; (size(find(Cons_Viol==0),1)/Np)*100];
    Constraint_violation=[Constraint_violation;best_Cons_Viol];
    
    
    
   if FEs>max_FEs
       break
   end
   if mod(FEs,10*D)==0
        Save_data.Best_Func_Val=[Save_data.Best_Func_Val; best_Func_Val];
        Save_data.Best_Cons_Viol=[Save_data.Best_Cons_Viol; best_Cons_Viol];
   end
  
end

if plot_flag == 1
    % Generation vs Objective Value Plotting %%%%%%%%%
    plot_Starting = 1;
    fig1 = figure(run);
    %         subplot(2,1,1);
    for i = plot_Starting:Gen-1
        if(Best_Cons_Viol(i) == 0)
            plot(i:i+1, Best_Func_Val(i:i+1), 'g.-');
            hold on;
        else
            plot(i:i+1, Best_Func_Val(i:i+1), 'r.-');
            hold on;
        end
    end
    title('Generation v/s Obj. Function Value of Best Solution');
    
    hold off;
    
    % Generation vs Feasibility(%) Plotting %%%%%%%%%
    fig2 = figure(1+run);
    for i = 1:Gen-1
        if(feasibility_Percentage(i) == 0)
            plot(i:i+1, feasibility_Percentage(i:i+1), 'r.-');
            hold on;
        else
            plot(i:i+1, feasibility_Percentage(i:i+1), 'g.-');
            hold on;
        end
    end
    title('Generation vs Feasibility(%)');
    hold on;
    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fig2 = figure(run+2);
    for i = 1:Gen-1
        if(Constraint_violation(i) == 0)
            plot(i:i+1, Constraint_violation(i:i+1), 'g.-');
            hold on;
        else
            plot(i:i+1, Constraint_violation(i:i+1), 'r.-');
            hold on;
        end
    end
    title('Generation vs Constrained of Violation');
    hold on;
    hold off;
end

end








