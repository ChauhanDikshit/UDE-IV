clear;clc;clear all

% User Defined Parameters %%%%%%%%%%%%%
for D=[30]
    dim=D;		%%%% Dimension %%%%
    Np = 100;
    frac_T = 0.25; % frac_T and T in the Unified_De file control ranking and modified ranking based parent selection
    frac_T1 = 0.25; % frac_T1 and T1 in the Unified_De file control number of sub-populations
    Lp = 25;
    max_FEs = 20000*D;
    H = 100;
    memory = 1;
    plot_flag = 0;
    
    para_adap = 1;
    sort_flag = 1;
    pc = 0.5;
    Max_gen = floor((max_FEs)/((frac_T*Np*3) + (1-frac_T)*Np));
    
    max_Run =25;
    
    global initial_flag
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Error_Vec = zeros(1,max_Run);
    
    all_Times = tic;
    
    for func_No = 1:28
        func_No;
        initial_flag = 0;
        
        rng(func_No);
        %rng('default');
        [x_Min, x_Max, numcon] = get_Bounds(func_No);
        
        Func1 = zeros(max_Run, 1);
        Cons1 = zeros(max_Run, 1);
        Func2 = zeros(max_Run, 1);
        Cons2 = zeros(max_Run, 1);
        Func3 = zeros(max_Run, 1);
        Cons3 = zeros(max_Run, 1);
        
        Solutions1 = zeros(max_Run, D);
        Solutions2 = zeros(max_Run, D);
        Solutions3 = zeros(max_Run, D);
        
        res_val1=[];
        Printing=1;
        
        for run = 1:max_Run
            parameters.run = run;
            parameters.D = D;
            parameters.func_No = func_No;
            parameters.Np = Np;
            parameters.x_Max = x_Max;
            parameters.x_Min = x_Min;
            parameters.max_Run = max_Run;
            parameters.frac_T = frac_T;
            parameters.frac_T1 = frac_T1;
            parameters.Lp = Lp;
            parameters.max_FEs = max_FEs;
            parameters.H = H;
            parameters.memory = memory;
            parameters.plot_flag = plot_flag;
            parameters.para_adap = para_adap;
            parameters.sort_flag = sort_flag;
            parameters.Max_gen = Max_gen;
            parameters.pc = pc;
            parameters.numcon = numcon;
            parameters.want_replacement=want_replacement;
            
            [Sol1, Sol2, Sol3,Save_data] = UDE_IV(parameters);
            
            % end
            Func1(run) = Sol1.func_Val;
            Cons1(run) = Sol1.cons_Viol;
            Solutions1(run, :) = Sol1.vector;
            
            Func2(run) = Sol2.func_Val;
            Cons2(run) = Sol2.cons_Viol;
            Solutions2(run, :) = Sol2.vector;
            
            Func3(run) = Sol3.func_Val;
            Cons3(run) = Sol3.cons_Viol;
            Solutions3(run, :) = Sol3.vector;
            
            str = sprintf('Func: %2d Run: %2d, Best Func Val: %f, Best Cons Val: %e', func_No, run, Func3(run), Cons3(run));
            disp(str);
            
            if Printing==1
                res_val= Save_data.Best_Func_Val;
                res_viol=Save_data.Best_Cons_Viol;
                res_val1=[res_val1 res_val,res_viol];
            end
        end
        
        Func = Func1;
        Cons = Cons1;
        Solutions = Solutions1;
        sorted_Id = sort1(Func, Cons);
        Result.Solutions = Solutions(sorted_Id, :);
        Result.Func_Val = Func(sorted_Id);
        Result.Cons_Val = Cons(sorted_Id);
        Result.Best = Result.Func_Val(1);
        Result.Median = Result.Func_Val(floor((max_Run+1)/2));
        Result.c = get_C(Result.Solutions(floor((max_Run+1)/2),:), func_No);
        Result.v_Bar = Result.Cons_Val(floor((max_Run+1)/2))/numcon;
        Result.Mean = mean(Result.Func_Val);
        Result.Worst = Result.Func_Val(end);
        Result.STD = std(Result.Func_Val);
        Result.SR = (size(find(Result.Cons_Val==0),1)/max_Run)*100;
        Result.vio_Bar = sum(Result.Cons_Val)/(max_Run*numcon);
        str = sprintf('2k_ResultFunc%d',func_No);
        %  	save(str, 'Result');
        
        Func = Func2;
        Cons = Cons2;
        Solutions = Solutions2;
        sorted_Id = sort1(Func, Cons);
        Result.Solutions = Solutions(sorted_Id, :);
        Result.Func_Val = Func(sorted_Id);
        Result.Cons_Val = Cons(sorted_Id);
        Result.Best = Result.Func_Val(1);
        Result.Median = Result.Func_Val(floor((max_Run+1)/2));
        Result.c = get_C(Result.Solutions(floor((max_Run+1)/2),:), func_No);
        Result.v_Bar = Result.Cons_Val(floor((max_Run+1)/2))/numcon;
        Result.Mean = mean(Result.Func_Val);
        Result.Worst = Result.Func_Val(end);
        Result.STD = std(Result.Func_Val);
        Result.SR = (size(find(Result.Cons_Val==0),1)/max_Run)*100;
        Result.vio_Bar = sum(Result.Cons_Val)/(max_Run*numcon);
        str = sprintf('10k_ResultFunc%d',func_No);
        %  	save(str, 'Result');
        
        Func = Func3;
        Cons = Cons3;
        Solutions = Solutions3;
        sorted_Id = sort1(Func, Cons);
        Result.Solutions = Solutions(sorted_Id, :);
        Result.Func_Val = Func(sorted_Id);
        Result.Cons_Val = Cons(sorted_Id);
        Result.Best = Result.Func_Val(1);
        Result.Median = Result.Func_Val(floor((max_Run+1)/2));
        Result.c = get_C(Result.Solutions(floor((max_Run+1)/2),:), func_No);
        Result.v_Bar = Result.Cons_Val(floor((max_Run+1)/2))/numcon;
        Result.Mean = mean(Result.Func_Val);
        Result.Worst = Result.Func_Val(end);
        Result.STD = std(Result.Func_Val);
        Result.SR = (size(find(Result.Cons_Val==0),1)/max_Run)*100;
        Result.vio_Bar = sum(Result.Cons_Val)/(max_Run*numcon);
        str = sprintf('20k_ResultFunc%d',func_No);
         % str=sprintf('20k_ResultFunc%d_d=%d',func_No,D);
        save(str, 'Result');
        
        str = sprintf('End of Function %d\n', func_No);
        disp(str);
        
        %% fitness values at different levels of the optimization process
        %%% required by the competition
        if Printing==1
            % lim=10*D:10*D:max_FEs;
            % res_to_print_val= res_val(:,lim);
            % res_to_print_viol= res_viol(:,lim);
            % res_to_print=[res_to_print_val;res_to_print_viol];
            % res_to_print=res_val1(:,lim);
            name1 = 'UDEIV_';
            name2 = num2str(func_No);
            name3 = '_';
            name4 = num2str(D);
            name5 = '.mat';
            f_name=strcat(name1,name2,name3,name4,name5);
            % res_to_print=res_to_print';
            save(f_name, 'res_val1', '-ascii','-append');
            name5 = '.dat';
            f_name=strcat(name1,name2,name3,name4,name5);
            % res_to_print=res_to_print';
            save(f_name, 'res_val1', '-ascii','-append');
            
        end
    end
    % getTable1;
    % getTable2;
     getTable3;    
end

