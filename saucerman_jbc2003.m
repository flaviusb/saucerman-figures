function output = saucerman_jbc2003(tspan)
% coupled signaling/EC for adult rat ventricular myocytes
%
% Reference:
% Jeffrey J. Saucerman, Laurence L. Brunton, Anushka P. Michailova, and Andrew D. McCulloch
% "Modeling beta-adrenergic control of cardiac myocyte contractility in silico", J. Biol. Chem., Vol 278: 47977-48003
%
% Copyright (2004) The Regents of the University of California
% All Rights Reserved
%
% Last modified: 4/19/2004
% Implemented by: Jeffrey Saucerman <jsaucer@ucsd.edu>
%
% Notes
% This version uses a reduced-order model of the L-type Ca channel: Michailova et al. (2004), submitted.
% For the original version, see the Berkeley Madonna code.

% Parameters
% ----- Signaling model parameters -------
% b-AR/Gs module
p(1) = 1;     % Ltotmax   [uM]
p(2) = 0.0132;  % sumb1AR   [uM]
p(3) = 3.83;    % Gstot     [uM]
p(4) = 0.285;   % Kl        [uM]
p(5) = 0.062;   % Kr        [uM]
p(6) = 33.0;    % Kc        [uM]
p(7) = 1.1e-3;  % k_barkp   [1/sec]
p(8) = 2.2e-3;  % k_barkm   [1/sec]
p(9) = 3.6e-3;  % k_pkap    [1/sec/uM]
p(10) = 2.2e-3; % k_pkam    [1/sec]
p(11) = 16.0;   % k_gact    [1/sec]
p(12) = 0.8;    % k_hyd     [1/sec]
p(13) = 1.21e3; % k_reassoc [1/sec/uM]
% cAMP module
p(14) = 49.7e-3;% AC_tot    [uM]
p(15) = 5.0e3;  % ATP       [uM]
p(16) = 38.9e-3;% PDEtot    [uM]
p(17) = 0.0;    % IBMXtot   [uM]
p(18) = 0.0;    % Fsktot    [uM] (10 uM when used)
p(19) = 0.2;    % k_ac_basal[1/sec]
p(20) = 8.5;    % k_ac_gsa  [1/sec]
p(21) = 7.3;    % k_ac_fsk  [1/sec]
p(22) = 5.0;    % k_pde     [1/sec]
p(23) = 1.03e3; % Km_basal  [uM]
p(24) = 315.0;  % Km_gsa    [uM]
p(25) = 860.0;  % Km_fsk    [uM]
p(26) = 1.3;    % Km_pde    [uM]
p(27) = 0.4;    % Kgsa      [uM]
p(28) = 44.0;   % Kfsk      [uM]
p(29) = 30.0;   % Ki_ibmx   [uM]
% PKA module
p(30) = 0.59;   % PKAItot   [uM]
p(31) = 0.025;  % PKAIItot  [uM]
p(32) = 0.18;   % PKItot    [uM]
p(33) = 9.14;   % Ka        [uM]
p(34) = 1.64;   % Kb        [uM]
p(35) = 4.375;  % Kd        [uM]
p(36) = 0.2e-3; % Ki_pki    [uM]
% PLB module
p(37) = 10;     % epsilon   [none]
p(38) = 106;    % PLBtot    [uM]
p(39) = 0.89;   % PP1tot    [uM]
p(40) = 0.3;    % Inhib1tot [uM]
p(41) = 54;     % k_pka_plb     [1/sec]
p(42) = 21;     % Km_pka_plb    [uM]
p(43) = 8.5;    % k_pp1_plb     [1/sec]
p(44) = 7.0;    % Km_pp1_plb    [uM]
p(45) = 60;     % k_pka_i1      [1/sec]
p(46) = 1.0;    % Km_pka_i1     [uM]
p(47) = 14.0;   % Vmax_pp2a_i1  [uM/sec]
p(48) = 1.0;    % Km_pp2a_i1    [uM]
p(49) = 1.0e-3; % Ki_inhib1     [uM]
% LCC module
p(50) = 25e-3;  % LCCtot        [uM]
p(51) = 25e-3;  % PP1lcctot     [uM]
p(52) = 25e-3;  % PP2Alcctot    [uM]
p(53) = 54;     % k_pka_lcc     [1/sec]
p(54) = 21;     % Km_pka_lcc    [uM]
p(55) = 8.52;   % k_pp1_lcc     [1/sec]
p(56) = 3;      % Km_pp1_lcc    [uM]
p(57) = 10.1;   % k_pp2a_lcc    [1/sec]
p(58) = 3;      % Km_pp2a_lcc   [uM]

% ---- EC Coupling model parameters ------
% universal parameters
p(59) = 20.8e-6;     % Vmyo  [uL]
p(60) = 9.88e-7;     % Vnsr  [uL]
p(61) = 9.3e-8;      % Vjsr  [uL]
p(62) = 1.534e-4;    % ACap  [cm^2] with C = 1 uF/cm^2
p(63) = 310;         % Temp  [K]
% extracellular concentrations     
p(64) = 140;     % Extracellular Na  [mM]
p(65) = 5.4;     % Extracellular K   [mM]
p(66) = 1.8;     % Extracellular Ca  [mM]
% current conductances
p(67) = 8.0;    % G_Na      [mS/uF]
p(68) = 0.35;   % G_to      [mS/uF] 
p(69) = 0.07;   % G_ss      [mS/uF]
p(70) = 0.24;   % G_kibar   [mS/uF] 
p(71) = 0.008;  % G_kp      [mS/uF]
% I_Ca parameters
p(72) = 300;    % f         [1/sec] 
p(73) = 2e3;    % g         [1/sec]
p(74) = 5187.5; % gammao    [1/sec/mM]
p(75) = 10;     % omega     [1/sec]
p(76) = 5.823e-9*3.0;  % pCa       [cm/sec]
p(77) = 1.078e-11*3.0;  % pK        [cm/sec]
p(78) = 3e5;    % Nlcc      [#/cell]
p(79) = -0.458; % I_Ca05    [uA/uF]
% pumps and background currents
p(80) = 1483;   % k_NaCa    [uA/uF]
p(81) = 87.5;   % Km_Na     [mM]
p(82) = 1.38;   % Km_Ca     [mM]
p(83) = 0.1;    % k_sat     [none]  
p(84) = 0.35;   % eta       [none]
p(85) = 1.1;  % ibarnak   [uA/uF]
p(86) = 10;     % Km_Nai    [mM]
p(87) = 1.5;    % Km_Ko     [mM]
p(88) = 1.15;   % ibarpca   [uA/uF]
p(89) = 0.5e-3; % Km_pca    [mM]
p(90) = 2.8e-3; % G_Cab     [uA/uF] 
p(91) = 1.18e-3;   % G_Nab     [uA/uF] 
p(92) = 0;      % Pns 
p(93) = 1.2e-3; % Km_ns     [mM]
% Calcium handling parameters
p(94) = 4.7;    % I_upbar   [mM/sec]
p(95) = 3e-4;   % Km_up     [mM]
p(96) = 15;     % nsrbar    [mM]
p(97) = 2e-3;   % tauon     [sec]
p(98) = 2e-3;   % tauoff    [sec]
p(99) = 60e3;   % gmaxrel   [mM/sec]
p(100) = 0.18e-3;% dcaith    [mM]
p(101) = 0.8e-3; % Km_rel    [mM]
p(102) = 8.75;   % CSQNth    [mM]
p(103) = 15;     % CSQNbar   [mM]
p(104) = 0.8;    % Km_csqn   [mM]
p(105) = 5.7e-4; % tau_tr    [sec]
p(106) = 0.07;   % TRPNbar   [mM]
p(107) = 0.05;   % CMDNbar   [mM]
p(108) = 0.07;   % INDObar   [mM]
p(109) = 0.5128e-3;  % Km_trpn   [mM]
p(110) = 2.38e-3;    % Km_cmdn   [mM]
p(111) = 8.44e-4;    % Km_indo   [mM]

% Initial conditions and mass matrix
% b-AR/Gsa module
%   1       2       3       4       5       6       7           8       9
%   L       R       Gs      b1ARtot b1ARd   b1ARp   Gsagtptot   Gsagdp  Gsbg
y10=[10;     0.01079;3.829;  0.01205;0.0;    1.154e-3;0.02505;   6.446e-4;0.02569];
m1 =[0,     0,      0,      1,      1,      1,      1,          1,      1];
% cAMP module and PKA module
%   10      11      12      13      14      15      
%   Gsa_gtp Fsk     AC      PDE     IBMX    cAMPtot 
y20=[0.02241;0;     0.04706;0.0389; 0;      0.8453];
m2 =[0,     0,      0,      0,      0,      1];
% PKA module
%   16      17      18
%   cAMP    PKACI   PKACII
y30=[0.2268;0.05868;8.278e-3];
m3 =[0,     0,      0];
% PLB 
%   19      20          21      22
%   PLBs    Inhib1ptot  Inhib1p PP1
y40=[4.105; 0.0526;     6.3e-5; 0.838];
m4 =[1,     1,          0,      0];
% LCC module
%   23          24
%   LCCap       LCCbp
y50=[5.103e-3;  5.841e-3];
m5 =[1,     1];
% Gating variables      
%   25      26      27      28      29      30      31      32      33
%   m       h       j       v       w       x       y       z       rto     sto     ssto    rss     sss   
y60=[1.4e-3;0.99;   0.99;   0.0;    0.0;    0.13;   0.96;   0.92;   1.4e-3; 1.0;    0.613;  198e-3; 0.43];
m6 =[1,     1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1,      1];
%   Intracellular concentrations/ Membrane voltage
%   38      39      40      41      42      43      44
%   Ca_nsr  Ca_jsr  Nai     Ki      Cai     Vm      trel
y70=[1.92;  1.92;   16;     145;    1.58e-4;-85.66; 0.9]; 
m7 =[1,     1,      1,      1,      1,      1,      1];

% Put everything together
y0  = [y10;y20;y30;y40;y50;y60;y70];    
M = diag([m1,m2,m3,m4,m5,m6,m7]); 

% Options
%tspan = [0;120];
options = odeset('Mass',M,'RelTol',1e-5,'MaxStep',5e-3,'Stats','on'); % set Reltol to 1e-6 for vclamp?

% Run simulation
[t,y] = ode15s(@f,tspan,y0,options,p);

yfinal = y(end,:)

Gsagtptot = y(:,7); cAMPtot = y(:,15); cAMPfree = y(:,16);
Nai = y(:,40); Ki = y(:,41); Cai = y(:,42); Vm = y(:,43);
Ca_nsr = y(:,38); Ca_jsr = y(:,39);
vlcc = y(:,28); wlcc = y(:,29); xlcc = y(:,30); ylcc = y(:,31); zlcc = y(:,32);

% global evars;
% evars = remove_repeats(evars);
% % Q_CaL = -1e3/10*trapz(evars(:,1),evars(:,2))   % Q_CaL [pC] should be 16.2+/-3
% % SRcontent = evars(end,3)
% te = evars(:,1); I_Ca = evars(:,2); I_rel = evars(:,3);
% clear global evars;

subplot(2,2,1);
plot(t,y(:,17));
title('Vm');

subplot(2,2,2);
plot(t,Cai);
title('Cai');
 
% subplot(2,2,3);
% plot(te,I_Ca);

output = [t,y];

function ydot = f(t,y,p)

ydot = zeros(size(y));
% -------- SIGNALING MODEL -----------

% b-AR module
LR = y(1)*y(2)/p(4);
LRG = LR*y(3)/p(5);
RG = y(2)*y(3)/p(6);
BARKDESENS = p(7)*(LR+LRG);
BARKRESENS = p(8)*y(5);
PKADESENS = p(9)*y(17)*y(4);  
PKARESENS = p(10)*y(6);
GACT = p(11)*(RG+LRG);
HYD = p(12)*y(7);
REASSOC = p(13)*y(8)*y(9);
ydot(1) = p(1)-LR-LRG-y(1);
ydot(2) = y(4)-LR-LRG-RG-y(2);
ydot(3) = p(3)-LRG-RG-y(3);
ydot(4) = (BARKRESENS-BARKDESENS)+(PKARESENS-PKADESENS);
ydot(5) = BARKDESENS-BARKRESENS;
ydot(6) = PKADESENS-PKARESENS;
ydot(7) = GACT-HYD;
ydot(8) = HYD-REASSOC;
ydot(9) = GACT-REASSOC;
% end b-AR module

% cAMP module
Gsa_gtp_AC = y(10)*y(12)/p(27);
Fsk_AC = y(11)*y(12)/p(28);
AC_ACT_BASAL = p(19)*y(12)*p(15)/(p(23)+p(15));	    
AC_ACT_GSA = p(20)*Gsa_gtp_AC*p(15)/(p(24)+p(15)); 
AC_ACT_FSK = p(21)*Fsk_AC*p(15)/(p(25)+p(15));	   
PDE_ACT = p(22)*y(13)*y(16)/(p(26)+y(16));	% FIXME: change 0.3*y(15) to cAMPfree   
PDE_IBMX = y(13)*y(14)/p(29);
ydot(10) = y(7)-Gsa_gtp_AC-y(10);
ydot(11) = p(18)-Fsk_AC-y(11);
ydot(12) = p(14)-Gsa_gtp_AC-y(12);  % note: assumes Fsk = 0.  Change Gsa_gtp_AC to Fsk_AC for Forskolin.
ydot(13) = p(16)-PDE_IBMX-y(13);
ydot(14) = p(17)-PDE_IBMX-y(14);
ydot(15) = AC_ACT_BASAL+AC_ACT_GSA+AC_ACT_FSK-PDE_ACT;
% end cAMP module

% PKA module
% 4/13/04 note: 
% I am currently not using an analytical solution for cAMPfree, but the
% below analytical solution does work.
% cAMPtemp = y(15)-(2*A2RC_I+2*A2R_I)-(2*A2RC_II+2*A2R_II);
% ydot(16) = cAMPtemp-sqrt(cAMPtemp^2-4*p(33)*(A2RC_I+A2RC_II))-y(16);
PKI = p(32)*p(36)/(p(36)+y(17)+y(18));
A2RC_I = (y(17)/p(35))*y(17)*(1+PKI/p(36));
A2R_I = y(17)*(1+PKI/p(36));
A2RC_II = (y(18)/p(35))*y(18)*(1+PKI/p(36));
A2R_II = y(18)*(1+PKI/p(36));
ARC_I = (p(33)/y(16))*A2RC_I;
ARC_II = (p(33)/y(16))*A2RC_II;
ydot(16) = y(15)-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-y(16);
PKAtemp = p(33)*p(34)/p(35)+p(33)*y(16)/p(35)+y(16)^2/p(35);
ydot(17) = 2*p(30)*y(16)^2-y(17)*(1+PKI/p(36))*(PKAtemp*y(17)+y(16)^2);
ydot(18) = 2*p(31)*y(16)^2-y(18)*(1+PKI/p(36))*(PKAtemp*y(18)+y(16)^2);
% end PKA module

% PLB module
PLB = p(38)-y(19);
PLB_PHOSPH = p(41)*y(17)*PLB/(p(42)+PLB);
PLB_DEPHOSPH = p(43)*y(22)*y(19)/(p(44)+y(19));
ydot(19) = PLB_PHOSPH-PLB_DEPHOSPH;
 
Inhib1 = p(40)-y(20);
Inhib1p_PP1 = y(21)*y(22)/p(49);
Inhib1_PHOSPH = p(45)*y(17)*Inhib1/(p(46)+Inhib1); 
Inhib1_DEPHOSPH = p(47)*y(20)/(p(48)+y(20));
ydot(20) = Inhib1_PHOSPH-Inhib1_DEPHOSPH;
ydot(21) = y(20)-Inhib1p_PP1-y(21);
ydot(22) = p(39)-Inhib1p_PP1-y(22);

fracPLBp = y(19)/p(38);
fracPLB = PLB/p(38);
fracPLBo = 0.9613; % adjust when changes are made to signaling model!
% end PLB module

% LCC module
LCCa = p(50)-y(23);
LCCa_PHOSPH = p(37)*p(53)*y(18)*LCCa/(p(54) + p(37)*LCCa);
LCCa_DEPHOSPH = p(37)*p(57)*p(52)*y(23)/(p(58)+p(37)*y(23));
ydot(23) = LCCa_PHOSPH - LCCa_DEPHOSPH;
fracLCCap = y(23)/p(50);
fracLCCapo = 0.2041;
 
LCCb = p(50)-y(24);
LCCb_PHOSPH = p(37)*p(53)*y(18)*LCCb/(p(54)+p(37)*LCCb);   
LCCb_DEPHOSPH = p(37)*p(55)*p(51)*y(24)/(p(56)+p(37)*y(24));
ydot(24) = LCCb_PHOSPH-LCCb_DEPHOSPH;
fracLCCbp = y(24)/p(50);
fracLCCbpo = 0.2336;
% end LCC module
% -------- END SIGNALING MODEL ---------
% -------- EC COUPLING MODEL -----------

% Constants
R = 8314;   % R     [J/kmol*K]
Frdy = 96485;  % Frdy     [C/mol]
FoRT = Frdy/R/p(63);
zna = 1;    % Na valence
zk = 1;     % K valence
zca = 2;    % Ca valence
% Nernst Potentials
ena = (1/FoRT/zna)*log(p(64)/y(40));       % should be 70.54 mV
ek = (1/FoRT/zk)*log(p(65)/y(41));		 % should be -87.94 mV
eca = (1/FoRT/zca)*log(p(66)/y(42));  % should be 120 mV
%eks = (1/FoRT)*logn((ko+prnak*nao)/(ki+prnak*nai))	% should be -77.54
ecl = -40.0;

% I_Na: Fast Na Current
am = 0.32*(y(43)+47.13)/(1-exp(-0.1*(y(43)+47.13)));
bm = 0.08*exp(-y(43)/11);
if y(43) >= -40
    ah = 0; aj = 0;
    bh = 1/(0.13*(1+exp(-(y(43)+10.66)/11.1)));
    bj = 0.3*exp(-2.535e-7*y(43))/(1+exp(-0.1*(y(43)+32)));
else
    ah = 0.135*exp((80+y(43))/-6.8);
    bh = 3.56*exp(0.079*y(43))+3.1e5*exp(0.35*y(43));
    aj = (-1.2714e5*exp(0.2444*y(43))-3.474e-5*exp(-0.04391*y(43)))*(y(43)+37.78)/(1+exp(0.311*(y(43)+79.23)));
    bj = 0.1212*exp(-0.01052*y(43))/(1+exp(-0.1378*(y(43)+40.14)));
end
ydot(25) = 1e3*(am*(1-y(25))-bm*y(25));
ydot(26) = 1e3*(ah*(1-y(26))-bh*y(26));
ydot(27) = 1e3*(aj*(1-y(27))-bj*y(27));
I_Na = p(67)*y(25)^3*y(26)*y(27)*(y(43)-ena);

% I_Ca: L-type Calcium Current
alcc = 400*exp( (y(43)+2)/10 );    % [1/sec]
blcc = 50*exp( -1*(y(43)+2)/13 );  % [1/sec]
flcc = p(72)*(0.375*fracLCCap/fracLCCapo+0.625); % PHOSPHOREGULATION
ylccinf = 1.0/(1+exp((y(43)+55.0)/7.5 )) +0.1/(1+exp((-y(43)+21.0)/6.0 ));
tauylcc = 0.02 + 0.3/( 1 +exp( (y(43)+30)/9.5 ) );    % [sec]			
gamma = p(74)*y(42); % [1/sec/mM]		
vgamma = gamma*((1-y(28))^4+2*y(28)*(1-y(28))^3+4*y(28)^2*(1-y(28))^2+8*y(28)^3*(1-y(28))+16*y(28)^4*(1-flcc/p(73)));
vomega = p(75)*((1-y(29))^4+1/2*y(29)*(1-y(29))^3+1/4*y(29)^2*(1-y(29))^2+1/8*y(29)^3*(1-y(29))+1/16*y(29)^4);

ydot(28) = alcc*(1-y(28))-blcc*y(28);
ydot(29) = 2*alcc*(1-y(29))-blcc/2*y(29);
ydot(30) = flcc*(1-y(30))-p(73)*y(30);
ydot(31) = (ylccinf - y(31))/tauylcc;
ydot(32) = vomega*(1-y(32))-vgamma*y(32);

if abs(y(43)) < 1e-6
    disp('Warning! Voltage near zero, could influence ibarca and ibark!');
end
ibarca = p(76)*4*(y(43)*Frdy*FoRT) * (1e-3*exp(2*y(43)*FoRT)-0.341*p(66)) /(exp(2*y(43)*FoRT)-1);
ibark = p(77)*(y(43)*Frdy*FoRT)*(y(41)*exp(y(43)*FoRT)-p(65)) /(exp(y(43)*FoRT)-1);
favail = 0.5*(0.4*fracLCCbp/fracLCCbpo+0.60);   % PHOSPHOREGULATION
I_Ca = ibarca*p(78)*favail*y(28)^4*y(30)*y(31)*y(32);
I_CaK = ibark/(1+I_Ca/p(79))*p(78)*favail*y(28)^4*y(30)*y(31)*y(32);
%I_Ca = 0; I_CaK = 0;
I_Catot = I_Ca+I_CaK;


% I_to: Transient Outward K Current
rtoss = 1/(1+exp((y(43)+10.6)/-11.42));
stoss = 1/(1+exp((y(43)+45.3)/6.8841));
taurto = 1.0/(45.16*exp(0.03577*(y(43)+50))+98.9*exp(-0.1*(y(43)+38))); % [sec]
tausto = 0.35*exp(-((y(43)+70)/15)^2)+0.035;	    % [sec]		 
taussto = 3.7*exp(-((y(43)+70)/30)^2)+0.035;        % [sec]			 
ydot(33) = (rtoss-y(33))/taurto;
ydot(34) = (stoss-y(34))/tausto;
ydot(35) = (stoss-y(35))/taussto;
I_to = p(68)*y(33)*(0.886*y(34)+0.114*y(35))*(y(43)-ek);   % [uA/uF]

% I_ss: Steady-state K Current
rssinf = 1/(1+exp(-(y(43)+11.5)/11.82));
taurss = 10/(45.16*exp(0.03577*(y(43)+50))+98.9*exp(-0.1*(y(43)+38)));
sssinf = 1/(1+exp((y(43)+87.5)/10.3));
tausss = 2.1;
ydot(36) = (rssinf-y(36))/taurss;
ydot(37) = (sssinf-y(37))/tausss;
I_ss = p(69)*y(36)*y(37)*(y(43)-ek);

% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(y(43)-ek-59.215)));
bki =(0.49124*exp(0.08032*(y(43)+5.476-ek)) + exp(0.06175*(y(43)-ek-594.31))) /(1 + exp(-0.5143*(y(43)-ek+4.753)));
kiss = aki/(aki+bki);
I_ki = p(70)*sqrt(p(65)/5.4)*kiss*(y(43)-ek) ;

% I_kp: Plateau K Current
kp = 1/(1+exp((7.488-y(43))/5.98));
I_kp = p(71)*kp*(y(43)-ek);

% I_ncx: Na/Ca Exchanger Current
s4 = exp(p(84)*y(43)*FoRT)*y(40)^3*p(66);
s5 = exp((p(84)-1)*y(43)*FoRT)*p(64)^3*y(42);
I_ncx = p(80)/(p(81)^3+p(64)^3) /(p(82)+p(66)) /(1+p(83)*exp((p(84)-1)*y(43)*FoRT)) *(s4-s5);

% I_nak: Na/K Pump Current
sigma = (exp(p(64)/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*y(43)*FoRT)+0.0365*sigma*exp(-y(43)*FoRT));
I_nak = p(85) *fnak /(1+(p(86)/y(40))^1.5) *(p(65)/(p(65)+p(87)));

% I_pca: Sarcolemmal Ca Pump Current
I_pca = p(88)*y(42)/(p(89)+y(42));
 
% I_cab: Ca Background Current
I_cab = p(90)*(y(43)-eca);
 
% I_nab: Na Background Current
I_nab = p(91)*(y(43)-ena);

% I_nsca: Nonspecific Ca-Activated Current: not used

% Total Membrane Currents
I_Na_tot = I_Na+I_nab+3*I_ncx+3*I_nak;          % [uA/uF]
I_K_tot = I_to+I_ss+I_ki+I_kp-2*I_nak+I_CaK;    % [uA/uF]
I_Ca_tot = I_Ca+I_cab+I_pca-2*I_ncx;            % [uA/uF]

% Calcium Induced Calcium Release (CICR)
trel = y(44)+2e-3;
ryron = 1-exp(-trel/p(97));  
ryroff = exp(-trel/p(98));
grel = p(99)/(1+exp((I_Ca_tot+5)/0.9));         % adjusted for rat  [1/sec] 
I_rel = grel*ryron*ryroff*(y(39)-y(42));        %   [mM/sec]

% Other SR fluxes and concentrations
Km_up = p(95)*(1+2*fracPLB)/(1+2*fracPLBo);     % PHOSPHOREGULATION
I_up = p(94)*y(42)^2/(Km_up^2+y(42)^2);         %   [mM/sec]
% Original: I_up = p(94)*y(42)^2/(p(95)^2+y(42)^2);         %   [mM/sec]
I_leak = p(94)*y(38)/p(96);                     %   [mM/sec]
I_tr = (y(38)-y(39))/p(105);                    %   [mM/sec]
Bjsr = 1/( 1+p(103)*p(104)/(p(104)+y(39))^2 );
ydot(38) = I_up-I_leak-I_tr*p(61)/p(60);        %   [mM/sec]
ydot(39) = Bjsr*(I_tr-I_rel);                   %   [mM/sec]
SRcontent = 1e3*((y(39)+y(39)/Bjsr)*p(61)/p(59)+y(38)*p(60)/p(59));    % [umol/L cytosol]

% Cytoplasmic Calcium Buffering
btrpn = p(106)*p(109)/(p(109)+y(42))^2;
bcmdn = p(107)*p(110)/(p(110)+y(42))^2;
bindo = p(108)*p(111)/(p(111)+y(42))^2;
Bmyo = 1/( 1+ bcmdn + btrpn + btrpn + bindo);

% Ion Concentrations and Membrane Potential
ydot(40) = -1e3*I_Na_tot*p(62)/(p(59)*zna*Frdy);          % [mM/sec] 
ydot(41) = -1e3*I_K_tot*p(62)/(p(59)*zk*Frdy);            % [mM/sec]
ydot(42) = -Bmyo*(1e3*I_Ca_tot*p(62)/(p(59)*zca*Frdy) ...
    +((I_up-I_leak)*p(60)/p(59))-(I_rel*p(61)/p(59)));    % [mM/sec]

% Simulation type
protocol = 'none';

switch lower(protocol)
    case {'none',''},
        I_app = 0;
    case 'pace',        % pace w/ current injection at rate 'rate'
		rate = 1;
		if mod(t+0.9,1/rate) <= 5e-3
            I_app = 10.0;
		else
            I_app = 0.0;
		end
    case 'vclamp',      
		V_hold = -40;
        V_test = -10;
		if (t > 59.1 & t < 59.5)
		    V_clamp = V_test;
		else
		    V_clamp = V_hold;
		end
		R_clamp = .02; 
		I_app = (V_clamp-y(43))/R_clamp;
end  

ydot(43) = -1e3*(I_Ca_tot+I_K_tot+I_Na_tot-I_app);

% CICR timing: y(44) tracks the last time when Vdot > 30 mV/msec
if (ydot(43) > 30e3) 
    ydot(44) = 1-1e4*y(44); 
else
    ydot(44) = 1;
end

% ----- END EC COUPLING MODEL ---------------

% Export intermediate variables (comment out for general usage)
% global evars;
% evars = [evars;t,I_Ca];%,SRcontent,I_nak];%,fracLCCap,fracLCCbp,fracPLB];