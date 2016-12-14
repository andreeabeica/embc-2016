%% run_stricker
%
% This script simulates the behaviour of the host-aware Stricker oscillator
% model.
%
% APS Darlington, a.p.s.darlington@warwick.ac.uk, December 2016
%

%% parameters

% host parameters
host.S0 = 1e4; host.sS = 0.5;
host.vT = 726; host.kT = 1000;
host.vE = 5800; host.kE = 1000;
host.nT = 300; host.nE = 300; host.nH = 300; host.nR = 7459;
host.wT = 4.14; host.wE = 4.14; host.wH = 948.93; host.wR = 930;
host.oT = 4.38; host.oE = 4.38; host.oH = 4.38; host.oR = 426.87; 
host.bX = 1; host.uX = 1;
host.kH = 152219; host.hH = 4;
host.dymX = 0.1; host.dypX = 0;
host.maxG = 1260; host.kG = 7;
host.M = 1e8; host.dyN = 0;

% circuit parameters
circ.w0 = 0.01;
circ.wA = 640; circ.wB = 320;
circ.oA = 4.38; circ.oB = 4.38;
circ.kA = 160; circ.kB = 160;
circ.hA = 2; circ.hB = 4;
circ.nA = 300; circ.nB = 300;
circ.bA = 1; circ.bB = 1;
circ.uA = 1; circ.uB = 1;
circ.dymA = 0.88; circ.dymB = 0.12;
circ.dypA = 0.29; circ.dypB = 0.06;

% inital conditions
iS0 = 0; ATP0 = 1000; 
mT0 = 0; cT0 = 0; pT0 = 0;
mE0 = 0; cE0 = 0; pE0 = 0;
mH0 = 0; cH0 = 0; pH0 = 0;
mR0 = 0; cR0 = 0; pR0 = 10;
mA0 = 0; cA0 = 0; pA0 = 10;
mB0 = 0; cB0 = 0; pB0 = 0;
N0 = 0;

Y0 = [iS0; ATP0; mT0; cT0; pT0; mE0; cE0; pE0; mH0; cH0; pH0; mR0; cR0; pR0; mA0; cA0; pA0; mB0; cB0; pB0; N0];

% time span
tmax = 1e4;

%% numerical integration of host-aware model

% simulate
[T,Y] = ode15s('stricker_host',[0,tmax],Y0,[],host,circ);

iS = Y(:,1); A = Y(:,2);
mT = Y(:,3); cT = Y(:,4); pT = Y(:,5);
mE = Y(:,6); cE = Y(:,7); pE = Y(:,8);
mH = Y(:,9); cH = Y(:,10); pH = Y(:,11);
mR = Y(:,12); cR = Y(:,13); pR = Y(:,14);
mA = Y(:,15); cA = Y(:,16); pA = Y(:,17);
mB = Y(:,18); cB = Y(:,19); pB = Y(:,20);
N = Y(:,21);

figure; plot(T,Y);

figure; plot(T(T>0.9*tmax),pA(T>0.9*tmax),T(T>0.9*tmax),pB(T>0.9*tmax)); legend('pA','pB');

