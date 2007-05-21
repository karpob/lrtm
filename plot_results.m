  figure(1)
  c=3.0e8;
  
  plot(100*(c./(f*1e9)),Tbeam_nadir,'b')
  hold on;
  plot(100*(c./(f*1e9)),Tbeam_limb,'r')
  xlabel('Wavelength (cm)')
  ylabel('T_B (^{\circ} K)')
  title('Tb')
  
  hold off;
  figure(2)
  plot(100*(c./(f*1e9)),R)
  xlabel('Wavelength (cm)')
  ylabel('Percent (%)')
  title('Limb Darkening')