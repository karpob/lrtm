clear all;
T=293;%:1:1000;
D=1e-6;
f=0:1:1000;
for i=1:length(f)
    [Np(i),Npp(i)]=refractivity_solution_cloud_93(T,D,1,f(i));
    tau(i)=3.336*Np(i);
    alpha(i)=0.182*Npp(i)*f(i);
end
for i=1:length(f)
    [Np2(i),Npp2(i)]=refractivity_solution_cloud_89(T,1e-6,1,f(i));
    tau2(i)=3.336*Np2(i);
    alpha2(i)=0.182*Npp2(i)*f(i);
end
for i=1:length(f)
    [Np3(i),Npp3(i)]=refractivity_solution_cloud_ulaby(T,1e-6,1,f(i));
    tau3(i)=3.336*Np3(i);
    alpha3(i)=0.182*Npp3(i)*f(i);
end

figure(1)
plot(f,tau)
hold on;
plot(f,tau2,'k')
plot(f,tau3, 'r')
legend('MPM 1993','MPM 1989 implemented','Ulaby \epsilon\prime \epsilon\prime\prime')
xlabel('Frequency (GHz)')
ylabel('Delay (ps/km)')
hold off;

figure(2);
plot(f,alpha)
hold on;
plot(f,alpha2,'k');
plot(f,alpha3,'r')
legend('MPM 1993', 'MPM 1989', 'Ulaby \epsilon\prime \epsilon\prime\prime')
xlabel('Frequency (GHz)')
ylabel('Absorption (dB/km)')
hold off;

figure(3)
plot(f,Np)
hold on;
plot(f,Np2,'k')
plot(f,Np3,'r')
legend('MPM 1993', 'MPM 1989', 'Ulaby \epsilon\prime \epsilon\prime\prime')
xlabel('Frequency (GHz)')
ylabel('N\prime')
hold off;

figure(4)
plot(f,Npp)
hold on;
plot(f,Npp2,'k')
plot(f,Npp3,'r')
legend('MPM 1993', 'MPM 1989', 'Ulaby \epsilon\prime \epsilon\prime\prime')
xlabel('Frequency (GHz)')
ylabel('N\prime\prime')
hold off;
