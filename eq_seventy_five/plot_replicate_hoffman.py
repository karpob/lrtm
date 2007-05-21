 #!/sw/bin/python2.3
from pylab import * 
matplotlib.rc('text',usetex= True)
matplotlib.use('PS')

figure(1) 
Tbeam_seventy_five=load('He_0_34_eq_Tbeam.dat')
plot(Tbeam_seventy_five,'+k')

#Tbeam_seventy_five_2=load('He_0_11_eq_Tbeam.dat')
#plot(Tbeam_seventy_five_2,'+b')

Tbeam_hoff_seventy_five=load('Tbeam_thesis_seventy_five.dat')
plot(Tbeam_hoff_seventy_five,'+r')

ylim(130,145) 

xticks(arange(15),('A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'))

legend(('2\,km Step 3.4\% He','Hoffman 2001'),loc='best',numpoints=1)

xlabel('Model Names')
ylabel('T ( $^{\circ}$K)')

title(r'$\theta_{zenith}$=75${^\circ}$ Brightness Temperatures')           

savefig('BT.ps',dpi=600)
 
figure(2)

residual=load('He_0_3_eq_residual.dat')
plot(residual,'+k')
 

#residual2=load('He_0_11_eq_residual.dat')
#plot(residual2,'+b')

ylim(0,8) 

xticks(arange(15),('A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'))

xlabel('Model Names')
ylabel('$\Delta$ T ($^{\circ}$K)')

#legend(('2\,km Step 3.4\% He','2\,km Step 11\% He'),loc='best',numpoints=1)
title('Difference vs. Hoffman 2001')           
savefig('res.ps',dpi=600)

figure(3)
title('Limb Darkening')

Tbeam1_nadir=load('Tbeam1_nadir.dat')
#Tbeam2_nadir=load('Tbeam2_nadir.dat')
Tbeam_hoff_nadir=load('Tbeam_thesis.dat')

plot(((Tbeam1_nadir-Tbeam_seventy_five)/Tbeam1_nadir)*100,'+k')
#plot(((Tbeam2_nadir-Tbeam_seventy_five_2)/Tbeam2_nadir)*100,'+b')
plot(((Tbeam_hoff_nadir-Tbeam_hoff_seventy_five)/Tbeam_hoff_nadir)*100,'+r')

ylim(0,10)
xlabel('Model Names')
ylabel('Limb Darkening (\%)')
xticks(arange(15),('A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'))
legend(('3.4\% He','Hoffman 2001'),loc='best',numpoints=1)
savefig('dark.ps',dpi=600)

figure(4)
title('Nadir $T_b$')
plot(Tbeam1_nadir,'+k')
#plot(Tbeam2_nadir,'+b')
plot(Tbeam_hoff_nadir,'+r')
legend(('3.4\% He','Hoffman 2001'),loc='best',numpoints=1)
xticks(arange(15),('A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'))
xlabel('Model Names')
ylabel('$T_b$')
savefig('nadir.ps',dpi=600)
