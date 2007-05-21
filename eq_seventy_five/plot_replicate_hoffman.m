
 figure(999)
 
 load seventy_five_eq_He_34
 dlmwrite('He_0_34_eq_Tbeam.dat',Tbeam_seventy_five,'\t')
 plot(Tbeam_seventy_five,'+k')
 hold on;
 
 load seventy_five_eq_He_11
 dlmwrite('He_0_11_eq_Tbeam.dat',Tbeam_seventy_five,'\t')
 plot(Tbeam_seventy_five,'+b')
 
 set(gca,'XTick',1:1:15)
 set(gca,'XTickLabel',Model_names)
 legend('2\,km Step 3.4\% He')
 title('Nadir Brightness Temperatures')           
 hold off
 
 figure(1000)
 load seventy_five_eq_He_34
 dlmwrite('He_0_34_eq_residual.dat',residual,'\t')
 plot(residual,'+b')
 hold on;
 load seventy_five_eq_He_11
 dlmwrite('He_0_11_eq_residual.dat',residual,'\t')
 plot(residual,'+k')
 
 set(gca,'XTick',1:1:15)
 set(gca,'XTickLabel',Model_names)
 legend('2\,km Step 3.4\% He')
 title('Nadir Brightness Differences Comparing against Hoffman')
 dlmwrite('Tbeam_thesis_seventy_five.dat',Tbeam_thesis_seventy_five,'\t')
 dlmwrite('Tbeam_thesis.dat',Tbeam_thesis,'\t')
 