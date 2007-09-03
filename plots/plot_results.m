  
  figure(1)
  c=3.0e8;
  load nominal.mat
  plot(100*(c./(f*1e9)),Tbeam_nadir,'b')
  hold on;
  load no_h2s.mat
  plot(100*(c./(f*1e9)),Tbeam_nadir,'r')
  load no_ph3.mat
  plot(100*(c./(f*1e9)),Tbeam_nadir,'m')
  load no_ph3_no_h2s.mat
  plot(100*(c./(f*1e9)),Tbeam_nadir,'k')
  
  xlabel('Wavelength (cm)')
  ylabel('T_B (^{\circ} K)')
  title('Tb')
  legend('Nominal','No H_2S','No PH_3', 'No PH_3 or H_2S')
  hold off;
  figure(2)
  
  load nominal.mat
  plot(100*(c./(f*1e9)),R,'b')
  hold on;
  load no_h2s.mat
  plot(100*(c./(f*1e9)),R,'r')
  load no_ph3.mat
  plot(100*(c./(f*1e9)),R,'m')
  load no_ph3_no_h2s.mat
  plot(100*(c./(f*1e9)),R,'k')
  hold on;
  legend('Nominal','No H_2S','No PH_3', 'No PH_3 or H_2S')
  
  xlabel('Wavelength (cm)')
  ylabel('Percent (%)')
  title('Limb Darkening')