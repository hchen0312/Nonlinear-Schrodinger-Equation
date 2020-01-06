% The algorithm includes Runge-Kutta method in its fourth order.
% Self-steepening and Raman scattering are involved, by giving proper
% estimation of Tr and taus

%---Specify input parameters
clear all; %
distance = 0.006; % in meter
loss_dB = -0; % in dB/cm
P_r=10^(loss_dB/10);
loss = -log(P_r)*100/2;
beta2 = -70*10^(-27); % in s2/m
beta3 = 3.6*10^(-40); % in s3/m, note that delta3=beta3/(6*beta2*T0)
beta4 = -11*10^(-56); % in s4/m
%beta5 = -0.3*10^(-69);
%beta6 = 4.1*10^(-84);
beta5=0;
beta6=0;

Aeff=0.96;
lambda=810*10^(-9);
gamma = 2*pi*23*10^(-20)/(lambda*Aeff*10^(-12)) % in 1/(w*m), this is the typical value in optical fiber
chirp0 = 0; % input pulse chirp (default value)
T0 = 100*10^(-15); % pulse width
V0=1/T0;
P0=3000; %in wat
N=sqrt(abs(gamma*P0*T0*T0/beta2/1.665/1.665))
Ld=abs(T0*T0/beta2)/1.665/1.665
Lnl=1/(gamma*P0)
L_fission=Ld/N
dleta_w=0.86*1/(T0/1.665)*gamma*P0*distance
f=3*10^(8)/lambda;
omega0=f*2*pi;
c=3*10^(8);
Tr=0*10^(-15);
taus=1/omega0;

%---set simulation parameters
nt = 256*2048/2; Tmax = 2*128*T0/2; % FFT points and window size
step_num = 4096/4; % No. of z steps
deltaz = distance/step_num; % step size in z
dtau = (2*Tmax)/nt; % step size in tau
%---Specify outer loop simulation paramters
nout=100;
dz=round(step_num/nout);
outer_count=1;
Map_t=zeros(nout-1,nt);
Map_f=zeros(nout-1,nt);
%---tau and omega arrays
tau = (-nt/2:nt/2-1)*dtau; % temporal grid
omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)]; % frequency grid
%---Input Field profile
uu = sqrt(P0)*exp(-1.39*(tau/T0).^2); %Gaussian shape
%---Plot input pulse shape and spectrum
temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi); % spectrum
figure; subplot(2,1,1);
plot (tau, abs(uu).^2, '--k'); hold on;
axis([-20*T0 20*T0 0 inf]);
xlabel('Time');
ylabel('Power');
title('Input and Output Pulse Shape and Spectrum');
subplot(2,1,2);
plot (c./(fftshift(omega)/(2*pi)+f)*10^(6), 10*log10(abs(temp).^2)-max(10*log10(abs(temp).^2)), '--k'); hold on;
axis([0.3,1.2 -40 5]);
%subplot(2,1,2);
%plot ((fftshift(omega)/(2*pi)+f)*10^(-12), 10*log10(abs(temp).^2)-max(10*log10(abs(temp).^2)), '--k'); hold on;
%axis([300,500 -30 5]);
xlabel('Frequency');
ylabel('Spectral Power');
%---Store dispersive phase shifts to speedup code
dispersion = exp(-0.5*loss*deltaz/2+0.5i*beta2*omega.^2*deltaz/2+1i*beta3/6*omega.^3*deltaz/2+1i*beta4/24*omega.^4*deltaz/2+1i*beta5/120*omega.^5*deltaz/2+1i*beta6/720*omega.^6*deltaz/2); % 1/2 phase factor
%hhz = 1i*gamma*deltaz/2; % 1/2 nonlinear phase factor
% ********* [ Beginning of MAIN Loop] ***********
% scheme: 4th order RK method

for n=1:step_num
ui=fft(ifft(uu).*dispersion);
%=======================================================================
diff1=1j*gamma*taus*gradient(abs(uu).^2.*uu)/dtau;
diff2=-gamma*Tr*uu.*gradient(abs(uu).^2)/dtau;
k1=1j*fft(ifft(gamma*abs(uu).*abs(uu).*uu+diff1+diff2).*dispersion);
%=======================================================================
uk2=ui+deltaz/2*k1;
diff3=1j*gamma*taus*gradient(abs(uk2).^2.*uk2)/dtau;
diff4=-gamma*Tr*uk2.*gradient(abs(uk2).^2)/dtau;
k2=1j*(gamma*abs(uk2).^2.*(uk2)+diff3+diff4);
%=======================================================================
uk3=ui+deltaz/2*k2;
diff5=1j*gamma*taus*gradient(abs(uk3).^2.*uk3)/dtau;
diff6=-gamma*Tr*uk3.*gradient(abs(uk3).^2)/dtau;
k3=1j*(gamma*abs(uk3).^2.*(uk3)+diff5+diff6);
%=======================================================================
uk4=fft(ifft(ui+deltaz*k3).*dispersion);
diff7=1j*gamma*taus*gradient(abs(uk4).^2.*uk4)/dtau;
diff8=-gamma*Tr*uk4.*gradient(abs(uk4).^2)/dtau;
k4=1j*(gamma*abs(uk4).^2.*(uk4)+diff7+diff8);
%=======================================================================
uu=fft(ifft(ui+k1*deltaz/6+k2*deltaz/3+k3*deltaz/3).*dispersion)+k4*deltaz/6;
temp = fftshift(ifft(uu)).* (nt*dtau)/sqrt(2*pi);
if n==outer_count*dz
    Map_t(outer_count,:)=abs(uu).^2;
    Map_f(outer_count,:)=10*log10(abs(temp).^2)-max(10*log10(abs(temp).^2));
    %Map_f(outer_count,:)=temp;
    outer_count=outer_count+1;
end

end
temp = fftshift(ifft(uu)).* (nt*dtau)/sqrt(2*pi); % Final spectrum
% *************** [ End of MAIN Loop ] **************
%----Plot output pulse shape and spectrum
subplot(2,1,1)
plot (tau, abs(uu).^2, '-k')
subplot(2,1,2)
plot(c./(fftshift(omega)/(2*pi)+f)*10^(6), 10*log10(abs(temp).^2)-max(10*log10(abs(temp).^2)), '-k')
%plot(c./(fftshift(omega)/(2*pi)+f)*10^(6), 10*log10(abs(temp).^2), '-k')
%subplot(2,1,2)
%plot((fftshift(omega)/(2*pi)+f)*10^(-12), 10*log10(abs(temp).^2)-max(10*log10(abs(temp).^2)), '-k')
output=10*log10(abs(temp).^2)-max(10*log10(abs(temp).^2));
output=output';
output_x=c./(fftshift(omega)/(2*pi)+f)*10^(6);
output_x=output_x';
%----Image Plot
hbar=1.055*10^(-34);
q=1.6*10^(-19);
shift_omega=fftshift(omega);
figure;
image([shift_omega(1)*hbar/q+1.24/(lambda*10^(6)),shift_omega(size(shift_omega))*hbar/q+1.24/(lambda*10^(6))],[distance,0],Map_f,'CDataMapping','scaled')
caxis([-40 0])
xlim([1,3])
figure;
image([-nt/2*dtau,(nt/2-1)*dtau],[0,distance],Map_t,'CDataMapping','scaled')
xlim([-20*T0,20*T0])
caxis([0,3000])