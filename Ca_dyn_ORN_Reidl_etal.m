function  [t,y,t2,k_Pact,tn,yn]=Ca_dyn_ORN_Reidl_etal(q);
% This matlab program integrates the four coupled nonlinear ordinary differential equations (12)-(15) in [1] for the special case of k_pCNGi1=k_pCNGi2.
% The parameters given below in the program code are taken from Tab.2 in [1].
% The output is both figures (directly produced by matlab) as well as a space-separated data file "concentrations.dat" containing the time and the concentrations of CNGo, Ca, CaM4 and CNGi.
% Run the program in MATLAB by typing Ca_dyn_ORN_Reidl_etal(q), where q determines the stimulus time pattern:
% q=2 increases the stimulus parameter k_pact twice for the durations dt and dt2 with a delay time (zeit2-zeit1) between the two pulses. Varying zeit2 leads to the results presented in Fig.3 in [1].
% q=3 increases k_pact once for 30s, leading to the data for Fig.4 in [1].
% q=4 again, increases k_pact twice, with varying time delay and/or different durations and stimulus strengths in the two pulses. Keeping the time delay fixed and varying the duration of the second pulse (same strengths as first) leads to data as presented in Fig.2.15 in [2].
% The data for Fig.2.10 in [2] can be produced using q=2 and only look at the first pulse, the data for Fig.2.11 in [2] with q=3 and varying stimulus strength.

% [1] J.Reidl, P.Borowski, A.Sensse, J.Starke, M.Zapotocky and M.Eiswirth, "Model of Calcium Oscillations Due to Negative Feedback in Olfactory Cilia", Biophys.J. 90, 1147 (2006).
% [2] P.Borowski, "Stochastic dynamics in olfactory signal transduction and development", Dissertation (2006), to be found at: http://nbn-resolving.de/urn:nbn:de:swb:14-1159519135136-22697

global t2
global k_Pact

t2=[];
k_Pact=[];
y0=[ 0 0 0 0 ];
tspan=[0 40000000000000];
options1=odeset('RelTol',10^(-8),'AbsTol',10^(-12),'OutputFcn','odeplot'); 
p=1; % getting initial conditions (rest conditions)
[t,y]=ode15s(@f2,tspan,y0,options1,p);
y0=[ y(end,1) y(end,2) y(end,3) y(end,4) ];

options1=odeset('RelTol',10^(-8),'AbsTol',10^(-12),'OutputFcn','odeplot','maxstep',.01); 
p=q; 
tspan=[-5 40]; % start and end time of simulation
datapoints=45*100*2;
if p==2
    tspan=[-2 14];
    datapoints=16*100*2;
end
if p==4
    tspan=[-2 4];
    datapoints=6*100*2;
end
[t,y]=ode15s(@f2,tspan,y0,options1,p);
    
figure(1)
plot(t,y(:,1));
set(1,'name','CNGo')
figure(2)
plot(t,y(:,2));
set(2,'name','Ca')
figure(3)
plot(t,y(:,3));
set(3,'name','CaM4')
figure(4)
plot(t,y(:,4));
set(4,'name','CNGi')

fid=fopen('concentrations.dat','wt');
for m=1:length(t)
    fprintf(fid,'%6.2e ',t(m));
    fprintf(fid,'%6.2e %6.2e %6.2e %6.2e\n',y(m,1),y(m,2),y(m,3),y(m,4));
end

function dydt=f2(t,y,p)
global t2
global k_Pact

% stimulus strengths:
k_pact=1.6e-5;  % resting activation of the system (no stimulus); p='plus'
reiz0=k_pact;   % no stimulation
reiz2=5.5;      % double pulse experiments
reiz3=5.5e-2;   % single, long pulse
reiz41=5.5;     % alternative two pulse exp., first pulse
reiz42=5.5;     % alternative two pulse exp., second pulse

% double pulse experiments (p=2):
% -------------------------------
% create two rectangle pulses with linear slopes
zeit1=0;    % start time first pulse
zeit2=4;    % start time second pulse
dt=.1;      % pulse duration first pulse
dt2=.1;     % pulse duration second pulse
ddt=.1;     % ratio of pulse duration used for slope at beginning and end of pulse

if t>zeit1 & p==2 & t<=zeit1+dt*ddt
k_pact=reiz0+(reiz2-reiz0)/(dt*ddt)*(t-zeit1);
end
if t>zeit1+dt*(1-ddt) & p==2 & t<=zeit1+dt
k_pact=reiz0+(reiz2-reiz0)/(dt*ddt)*(-t+zeit1+dt);
end
if t>zeit2 & p==2 & t<=zeit2+dt2*ddt
k_pact=reiz0+(reiz2-reiz0)/(dt2*ddt)*(t-zeit2);
end
if t>zeit2+dt2*(1-ddt) & p==2 & t<=zeit2+dt2
k_pact=reiz0+(reiz2-reiz0)/(dt2*ddt)*(-t+zeit2+dt2);
end
if t>zeit1+dt*ddt & p==2 & t<=zeit1+dt*(1-ddt)
k_pact=reiz2;
end
if t>zeit2+dt2*ddt & p==2 & t<=zeit2+dt2*(1-ddt)
k_pact=reiz2;
end

% single long pulse (p=3):
% ------------------------
zeit3=0;                    % start time for long pulse
if t>zeit3 & p==3 & t<30    % pulse duration
k_pact=reiz3;
end

% two pulses, different strengths and durations (p=4):
% ----------------------------------------------------
pulseduration1=.2;
pulseduration2=.1;
ptime1=0;
ptime2=2.5;
if t>ptime1 & p==4 & t<ptime1+pulseduration1
    k_pact=reiz41;
end
if t>ptime2 & p==4 & t<ptime2+pulseduration2
    k_pact=reiz42;
end

% parameters (cf. Tab.2 in [1])
sigma=5E-7;         % volume/surface ratio of the cilium
CNG_tot=1.3E-13;    % surface concentration of CNG channels
CaM_tot=2E-5;       % concentration of calmodulin
k_mCNGo= 1e-2;      % m='minus'
k_mCaM4=2.5;
k_pCaM4=1.1e9;      % p='plus'
k_pCNGi=2.1e6;
k_mCNGi=3.4e-1;
i_Ca=2e4;           % calcium current through a single CNG channel
k_Ca=1e-10;         % calcium extrusion
K_Ca=1.2e-7;        % calcium extrusion

% variables
CNGo=y(1);          % surface concentration of open channels
Ca=y(2);            % calcium concentration
CaM4=y(3);          % concentration of calmodulin-calcium complex
CNGi=y(4);          % surface concentration of inhibited channels

% writing out the simulus time course (optional)
if (p==2 | p==3 | p==4) 
    t2=[t2 t];
    k_Pact=[k_Pact k_pact];
end

% differential equations (k_pCNGi1=k_pCNGi2=k_pCNGi)

dCNGo=k_pact*(CNG_tot-CNGo-CNGi)-k_pCNGi*CaM4*CNGo-k_mCNGo*CNGo;

dCaM4=k_pCaM4*Ca^2*(CaM_tot-CaM4-CNGi/sigma)-k_mCaM4*CaM4-k_pCNGi/sigma*CaM4*(CNG_tot-CNGi)+k_mCNGi/sigma*CNGi;

dCa=CNGo/sigma*i_Ca-k_Ca/sigma*Ca/(Ca+K_Ca)-4*(k_pCaM4*Ca^2*(CaM_tot-CaM4-CNGi/sigma)-k_mCaM4*CaM4);

dCNGi=-k_mCNGi*CNGi+k_pCNGi*CaM4*(CNG_tot-CNGi);

% rewriting the differentials in matlab suitable form
dydt(1,1)=dCNGo;
dydt(2,1)=dCa;
dydt(3,1)=dCaM4;
dydt(4,1)=dCNGi;
