clear
 figure(999)
 
 load He_0_3_frozen_eq
 dlmwrite('He_0_3_frozen_eq_Tbeam.dat',Tbeam_nadir,'\t')
 plot(Tbeam_nadir,'+b')
 hold on;
 
 load He_0_34_frozen_eq
 dlmwrite('He_0_34_frozen_eq_Tbeam.dat',Tbeam_nadir,'\t')
 plot(Tbeam_nadir,'+k')
 
 load He_0_9_frozen_eq
 dlmwrite('He_0_9_frozen_eq_Tbeam.dat',Tbeam_nadir,'\t')
 plot(Tbeam_nadir,'+g')
 
 load He_0_11_frozen_eq
 dlmwrite('He_0_11_frozen_eq_Tbeam.dat',Tbeam_nadir,'\t')
 plot (Tbeam_nadir,'+r')
 
 load He_0_3_frozen_eq_auto
 dlmwrite('He_0_3_frozen_eq_auto_Tbeam.dat',Tbeam_nadir,'\t')
 plot (Tbeam_nadir, 'ob')

 load He_0_34_frozen_eq_auto
 dlmwrite('He_0_34_frozen_eq_auto_Tbeam.dat',Tbeam_nadir,'\t')
 plot (Tbeam_nadir,'ok')
 
  
 load He_0_9_frozen_eq_auto
 dlmwrite('He_0_9_frozen_eq_auto_Tbeam.dat',Tbeam_nadir,'\t')
 plot (Tbeam_nadir, 'og')
 
 load He_0_11_frozen_eq_auto
 dlmwrite('He_0_11_frozen_eq_auto_Tbeam.dat',Tbeam_nadir,'\t')
 plot (Tbeam_nadir, 'or')
 
 set(gca,'XTick',1:1:15)
 set(gca,'XTickLabel',Model_names)
 legend('2\,km Step 3% He','2\,km Step 3.4% He','2\,km Step 9% He','2\,km Step 11% He', ...
        'AutoStep 3% He','AutoStep 3.4% He','AutoStep 9% He','AutoStep 11% He')
 title('Nadir Brightness Temperatures')           
 hold off
 
 figure(1000)
 load He_0_3_frozen_eq
 dlmwrite('He_0_3_frozen_eq_residual.dat',residual,'\t')
 plot(residual,'+b')
 
 hold on;
 
 load He_0_34_frozen_eq
 dlmwrite('He_0_34_frozen_eq_residual.dat',residual,'\t')
 plot(residual,'+k')
 
 load He_0_9_frozen_eq
 dlmwrite('He_0_9_frozen_eq_residual.dat',residual,'\t')
 plot(residual,'+g')
 
 load He_0_11_frozen_eq
 dlmwrite('He_0_11_frozen_eq_residual.dat',residual,'\t')
 plot (residual,'+r')

 load He_0_3_frozen_eq_auto
 dlmwrite('He_0_3_frozen_eq_auto_residual.dat',residual,'\t')
 plot (residual, 'ob')
 
 load He_0_34_frozen_eq_auto
 dlmwrite('He_0_34_frozen_eq_auto_residual.dat',residual,'\t')
 plot (residual,'ok')
 
 load He_0_9_frozen_eq_auto
 dlmwrite('He_0_9_frozen_eq_auto_residual.dat',residual,'\t')
 plot (residual, 'og')
 
 load He_0_11_frozen_eq_auto
 dlmwrite('He_0_11_frozen_eq_auto_residual.dat',residual,'\t')
 plot (residual, 'or')
 
 set(gca,'XTick',1:1:15)
 set(gca,'XTickLabel',Model_names)
 legend('2\,km Step 3% He','2\,km Step 3.4% He','2\,km Step 9% He','2\,km Step 11% He', ...
        'AutoStep 3% He','AutoStep 3.4% He','AutoStep 9% He','AutoStep 11% He')
 title('Nadir Brightness Differences Comparing against Hoffman')           
 hold off
 
 
 
 
 
% 
% figure(1000)
% plot(Tbeam_seventy_five-Tbeam_thesis_seventy_five,'+k')
% hold on;
% plot(Tbeam_nadir-Tbeam_thesis,'+b')
% set(gca,'XTickLabel',Model_names)
% set(gca,'XTick',1:1:15)
% set(gca,'XTickLabel',Model_names)
% legend('Error in Limb values ^{\circ} K','Error in nadir values ^{\circ} K')
% title('Error compared with Hoffman')
% %savefig('error','pdf','-cmyk','-r600')
% hold off
% 
% figure(1001)
% plot(Tbeam_thesis-Tbeam_thesis_seventy_five,'+r')
% hold on
% plot(Tbeam_nadir-Tbeam_seventy_five,'+g')
% set(gca,'XTick',1:1:15)
% set(gca,'XTickLabel',Model_names)
% legend('Hoffman Values','Run values')
% title('Compare Limb darkening')
% %savefig('darkening','pdf','-cmyk','-r600')
% hold off
% %mv *.pdf plots
% /*P T dr XH2 XHe XH2S XNH3 XH2O XCH4 XPH3  clouds DNH4SH DH2S DNH3 DH2O DCH4
% DPH3 DSOL g mu refr_w/o refr_ w/  */
color_funtime={'b','g','r','k','y','--b','--g','--r','--k','--y','+b','+g','+r','+k','+y'}
figure(1002)
for j=1:1
        plot(XPH3(1:296,1),tcme(1:296,1,1),color_funtime{j})
        hold on;
end
hold off