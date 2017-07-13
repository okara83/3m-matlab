function pppos = probepos()

% 9 probes at the equator corresponding to (2*4+1) = 10 points minus g_0^0 
%for m_max=4, then there is [4,7,9,7,4] probes respectively symetric to the
%equator.
% Colatitude correspond to the 5 Gauss's points to solve l_max=4.

%reordered 101214 for the actual probe ordering!

phi_HLN=[22.59;	112.72;	202.63;	292.36];
phi_MLN=[25.7;77.11;128.6;180;231.4;282.9;334.29];
phi_EQ=[0;	39.99;	80.01;	120.06;	159.63;	200.01;	240.1;	280.03;	320];
phi_MLS=[25.7;	77.11;	128.6;	180;	231.4;	282.9;	334.31];
phi_HLS=[45.01;	135.01;	224.94;	315];
phi=[phi_HLS;phi_MLS;phi_EQ;phi_MLN;phi_HLN];           %changed by Onur - was upside down
theta_HLN=[53.01;	51.54;	51.06;	53.51];					
theta_MLN=[28.02;	27.47;	26.84;	26.55;	25.96;	25.14;	28.31];
theta_EQ=[-2.58;	-2.95;	-3.51;	-4.11;	-0.93;	-1.31;	-1.72;	-1.72;	-2.35];
theta_MLS=[-32.78;	-33.13;	-33.80;	-34.39;	-35.21;	-35.71;	-32.47];	
theta_HLS=[-65.07;	-64.98;	-64.98;	-65.97];		
theta=[theta_HLS;theta_MLS;theta_EQ;theta_MLN;theta_HLN]-90;

r(1,1:31)=1.052191781; %in units of "core mantle boundary" radius of 1.46m.  r (m) = 1.536m


%% %%%%%%%%%%%%%%%%%%% Uncomment to plot the probes position %%%%%%%%%%%%%%
% figure('name','probe positions')
% plot(phi,abs(theta),'or'); 
% xlabel('longitude (degree)'); 
% ylabel('Colatitude (degree)'); 
% title('Positions of the Hall probes');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert in radian
theta=abs(theta)*pi/180;
phi=phi*pi/180;

pppos=[r',theta,phi]
