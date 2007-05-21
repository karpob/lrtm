 #!/sw/bin/python2.3
from pylab import * 
matplotlib.rc('text',usetex= True)
matplotlib.use('PS')
Tbeam_nadir=load('He_0_3_frozen_eq_Tbeam.dat')
plot(Tbeam_nadir,'+b')
#hold on;
 
Tbeam_nadir=load('He_0_34_frozen_eq_Tbeam.dat')
plot(Tbeam_nadir,'+k')
 
Tbeam_nadir=load('He_0_9_frozen_eq_Tbeam.dat')
plot(Tbeam_nadir,'+g')
 
 
Tbeam_nadir=load('He_0_11_frozen_eq_Tbeam.dat')
plot (Tbeam_nadir,'+r')
 
 
Tbeam_nadir=load('He_0_3_frozen_eq_auto_Tbeam.dat')
plot (Tbeam_nadir, 'xb')

Tbeam_nadir=load('He_0_34_frozen_eq_auto_Tbeam.dat')
plot (Tbeam_nadir,'xk')
 
Tbeam_nadir=load('He_0_9_frozen_eq_auto_Tbeam.dat')
plot (Tbeam_nadir, 'xg')
 
Tbeam_nadir=load('He_0_11_frozen_eq_auto_Tbeam.dat')
plot (Tbeam_nadir, 'xr')
 
#set(gca,'XTick',:)
xticks(arange(15),('A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'))
legend(('2\,km Step 3\% He','2\,km Step 3.4\% He','2\,km Step 9\% He','2\,km Step 11\% He','AutoStep 3\% He','AutoStep 3.4\% He','AutoStep 9\% He','AutoStep 11\% He'),loc='best',numpoints=1)
xlabel('Model Names')
ylim(140,150)
ylabel('T ( $^{\circ}$K)')
title('Nadir Brightness Temperatures')           
#hold off
savefig('BT.ps',dpi=600)
 
figure(2)
#load He_0_3_frozen_eq
residual=load('He_0_3_frozen_eq_residual.dat')
plot(residual,'+b')
 
#hold on;

#load He_0_34_frozen_eq
residual=load('He_0_34_frozen_eq_residual.dat')
plot(residual,'+k')
 
#load He_0_9_frozen_eq
residual=load('He_0_9_frozen_eq_residual.dat')
plot(residual,'+g')
 
#load He_0_11_frozen_eq
residual=load('He_0_11_frozen_eq_residual.dat')
plot (residual,'+r')

#load He_0_3_frozen_eq_auto
residual=load('He_0_3_frozen_eq_auto_residual.dat')
plot (residual, 'xb')
 
#load He_0_34_frozen_eq_auto
residual=load('He_0_34_frozen_eq_auto_residual.dat')
plot (residual,'xk')
 
#load He_0_9_frozen_eq_auto
residual=load('He_0_9_frozen_eq_auto_residual.dat')
plot (residual, 'xg')
 
#load He_0_11_frozen_eq_auto
residual=load('He_0_11_frozen_eq_auto_residual.dat')
plot (residual, 'xr')
 
#set(gca,'XTick',1:1:15)
xticks(arange(15),('A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'))
xlabel('Model Names')
ylabel('$\Delta$ T ($^{\circ}$K)')
ylim(0,8)
#set(gca,'XTickLabel',Model_names)
legend(('2\,km Step 3\% He','2\,km Step 3.4\% He','2\,km Step 9\% He','2\,km Step 11\% He','AutoStep 3\% He','AutoStep 3.4\% He','AutoStep 9\% He','AutoStep 11\% He'),loc='best',numpoints=1)
title('Difference vs. Hoffman 2001')           
savefig('res.ps',dpi=600)
