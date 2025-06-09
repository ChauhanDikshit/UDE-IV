cd 'C:\Users\diksh_06\OneDrive - National University of Singapore\Dikshit Chauhan NUS\Matlab code of UDE-III\Stagnated Plots\D=30 with and without replacement'
for i=1:28
fig1=figure(i);func_No=i;
legend('With Replacement','Fontname','Times New Roman','FontSize',12,'FontWeight','bold')
 legend('Without Replacement','With Replacement','Fontname','Times New Roman','FontSize',12,'FontWeight','bold')
 saveas(fig1, strcat(strcat('Centroid population','_Func_num_',num2str(func_No), '_D_',num2str(D),'_run_',num2str(run)),'.png'));
 saveas(fig1, strcat(strcat('Centroid population','_Func_num_',num2str(func_No), '_D_',num2str(D),'_run_',num2str(run)),'.fig'));
end