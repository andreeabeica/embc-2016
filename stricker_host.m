function dY = stricker_host(T,Y,~,host, circ);

%% stricker_host
%
% This function contains the ODEs describing the Stricker oscillator
% architecture embedded within the Weisse et al host model.
%
% APS Darlington, a.p.s.darlington@warwick.ac.uk, December 2016
% 

%% define species
iS = Y(1); ATP = Y(2);
mT = Y(3); cT = Y(4); pT = Y(5);
mE = Y(6); cE = Y(7); pE = Y(8);
mH = Y(9); cH = Y(10); pH = Y(11);
mR = Y(12); cR = Y(13); pR = Y(14);
mA = Y(15); cA = Y(16); pA = Y(17); % activator
mB = Y(18); cB = Y(19); pB = Y(20); % repressor
N = Y(21);

%% host parameters
eS = host.S0;
vT = host.vT; kT = host.kT;
vE = host.vE; kE = host.kE;
sS = host.sS;
nT = host.nT; nE = host.nE; nH = host.nH; nR = host.nR;
wT = host.wT; wE = host.wE; wH = host.wH; wR = host.wR;
oT = host.oT; oE = host.oE; oH = host.oH; oR = host.oR;
kH = host.kH; hH = host.hH;
bX = host.bX; uX = host.uX;
M = host.M;
dymX = host.dymX; dypX = host.dypX;
maxG = host.maxG;
kG = host.kG;
dyN = host.dyN;

%% circuit parameters

% energy-independent leakiness
w0 = circ.w0;

% circuit parameters
wA = circ.wA; wB = circ.wB;
oA = circ.oA; oB = circ.oB;
kA = circ.kA; kB = circ.kB;
hA = circ.hA; hB = circ.hB;
nA = circ.nA; nB = circ.nB;
bA = circ.bA; bB = circ.bB;
uA = circ.uA; uB = circ.uB;
dymA = circ.dymA; dymB = circ.dymB;
dypA = circ.dypA; dypB = circ.dypB;

%% growth rate
gamma = (maxG*ATP)/(kG+ATP);
translatingribosomes = cT + cE + cH + cR + cA + cB;
lambda = (1/M)*gamma*translatingribosomes;

%% rates
g2mT = (wT*ATP)/(oT + ATP);
g2mE = (wE*ATP)/(oE + ATP);
g2mH = ((wH*ATP)/(oH + ATP))*(1/(1+(pH/kH)^hH));
g2mR = (wR*ATP)/(oR + ATP);
g2mA = w0 + ((wA*ATP)/(oA + ATP))*(((pA/kA)^hA)/(1+((pA/kA)^hA)+((pB/kB)^hB)+((pA/kA)^hA)*((pB/kB)^hB)));
g2mB = w0 + ((wB*ATP)/(oB + ATP))*(((pA/kA)^hA)/(1+((pA/kA)^hA)+((pB/kB)^hB)+((pA/kA)^hA)*((pB/kB)^hB)));

m2pT = (cT*gamma)/nT;
m2pE = (cE*gamma)/nE;
m2pH = (cH*gamma)/nH;
m2pR = (cR*gamma)/nR;
m2pA = (cA*gamma)/nA;
m2pB = (cB*gamma)/nB;

%% host ODEs
% metabolism
diS = (pT*(vT*eS)/(kT+eS)) - (pE*(vE*iS)/(kE+iS)) - lambda*iS;
dATP = (sS*pE*(vE*iS)/(kE+iS)) - lambda*ATP...
    - nR*m2pR - nT*m2pT - nE*m2pE - nH*m2pH...
    - nA*m2pA - nB*m2pB;

% transport
dmT = g2mT - (lambda+dymX)*mT + m2pT - bX*pR*mT + uX*cT;
dcT = -lambda*cT + bX*pR*mT - uX*cT - m2pT;
dpT = m2pT - (lambda+dypX)*pT;

% enzymes
dmE = g2mE - (lambda+dymX)*mE + m2pE - bX*pR*mE + uX*cE;
dcE = -lambda*cE + bX*pR*mE - uX*cE - m2pE;
dpE = m2pE - (lambda+dypX)*pE;

% host proteins
dmH = g2mH - (lambda+dymX)*mH + m2pH - bX*pR*mH + uX*cH;
dcH = -lambda*cH + bX*pR*mH - uX*cH - m2pH;
dpH = m2pH - (lambda+dypX)*pH;

% host ribosomes
dmR = g2mR - (lambda+dymX)*mR + m2pR - bX*pR*mR + uX*cR;
dcR = -lambda*cR + bX*pR*mR - uX*cR - m2pR;
dpR = m2pR - (lambda+dypX)*pR...
    + m2pE - bX*pR*mE + uX*cE...
    + m2pT - bX*pR*mT + uX*cT...
    + m2pH - bX*pR*mH + uX*cH...
    + m2pR - bX*pR*mR + uX*cR...
    + m2pA - bA*pR*mA + uA*cA...
    + m2pB - bB*pR*mB + uB*cB;

%% circuit ODEs
dmA = g2mA - (lambda+dymA)*mA + m2pA - bA*pR*mA + uA*cA;
dcA = -lambda*cA + bA*pR*mA - uA*cA - m2pA;
dpA = m2pA - (lambda+dypA)*pA;

dmB = g2mB - (lambda+dymB)*mB + m2pB - bB*pR*mB + uB*cB;
dcB = -lambda*cB + bB*pR*mB - uB*cB - m2pB;
dpB = m2pB - (lambda+dypB)*pB;

%% population
dN = (lambda-dyN)*N;

%% return dY
dY = [diS; dATP; dmT; dcT; dpT; dmE; dcE; dpE; dmH; dcH; dpH; dmR; dcR; dpR; dmA; dcA; dpA; dmB; dcB; dpB; dN];

end
