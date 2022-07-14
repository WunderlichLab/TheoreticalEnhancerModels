% function [model] = MEC_2_ZC_2_c_enhancer_13_syms(f0_user)
function [model] = MEC_2_ZC_2_c_enhancer_13_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms mu_1 mu_2 mu_3 mu_4 mu_5 mu_6 mu_7 C_1_1 C_1_2 C_1_3 C_1_4 C_1_5 C_1_6 C_1_7 C_2_2 C_2_3 C_2_4 C_2_5 C_2_6 C_2_7 C_3_3 C_3_4 C_3_5 C_3_6 C_3_7 C_4_4 C_4_5 C_4_6 C_4_7 C_5_5 C_5_6 C_5_7 C_6_6 C_6_7 C_7_7

x = [
mu_1, mu_2, mu_3, mu_4, mu_5, mu_6, mu_7, C_1_1, C_1_2, C_1_3, C_1_4, C_1_5, C_1_6, C_1_7, C_2_2, C_2_3, C_2_4, C_2_5, C_2_6, C_2_7, C_3_3, C_3_4, C_3_5, C_3_6, C_3_7, C_4_4, C_4_5, C_4_6, C_4_7, C_5_5, C_5_6, C_5_7, C_6_6, C_6_7, C_7_7 ...
];

% PARAMETERS

syms kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 

% KAPPA (constant parameters)

syms Omega indmu1 indmu2 indmu3 indmu4 indmu5 indmu6 indmu7 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 indC11 indC12 indC13 indC14 indC15 indC16 indC17 indC18 indC19 indC20 indC21 indC22 indC23 indC24 indC25 indC26 indC27 indC28 kmu01 kmu02 kmu03 kmu04 kmu05 kmu06 kmu07 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 kC011 kC012 kC013 kC014 kC015 kC016 kC017 kC018 kC019 kC020 kC021 kC022 kC023 kC024 kC025 kC026 kC027 kC028 

syms t

p = [kon1,koff1,r1,kon2,koff2,r2,myalpha,beta1,betam1,beta2,betam2];

k = [Omega,indmu1,indmu2,indmu3,indmu4,indmu5,indmu6,indmu7,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,indC11,indC12,indC13,indC14,indC15,indC16,indC17,indC18,indC19,indC20,indC21,indC22,indC23,indC24,indC25,indC26,indC27,indC28,kmu01,kmu02,kmu03,kmu04,kmu05,kmu06,kmu07,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010,kC011,kC012,kC013,kC014,kC015,kC016,kC017,kC018,kC019,kC020,kC021,kC022,kC023,kC024,kC025,kC026,kC027,kC028];

if nargin > 0
   f0_user = varargin{1};
   if ~isnumeric(f0_user)
      p_user = setdiff(symvar(f0_user),p);
      % ADDITIONAL PARAMETERS IN INITIAL CONDITIONS
      p = [p,p_user];
   end
	fmu01 = f0_user(1); 
	fmu02 = f0_user(2); 
	fmu03 = f0_user(3); 
	fmu04 = f0_user(4); 
	fmu05 = f0_user(5); 
	fmu06 = f0_user(6); 
	fmu07 = f0_user(7); 
	fC01 = f0_user(8); 
	fC02 = f0_user(9); 
	fC03 = f0_user(10); 
	fC04 = f0_user(11); 
	fC05 = f0_user(12); 
	fC06 = f0_user(13); 
	fC07 = f0_user(14); 
	fC08 = f0_user(15); 
	fC09 = f0_user(16); 
	fC010 = f0_user(17); 
	fC011 = f0_user(18); 
	fC012 = f0_user(19); 
	fC013 = f0_user(20); 
	fC014 = f0_user(21); 
	fC015 = f0_user(22); 
	fC016 = f0_user(23); 
	fC017 = f0_user(24); 
	fC018 = f0_user(25); 
	fC019 = f0_user(26); 
	fC020 = f0_user(27); 
	fC021 = f0_user(28); 
	fC022 = f0_user(29); 
	fC023 = f0_user(30); 
	fC024 = f0_user(31); 
	fC025 = f0_user(32); 
	fC026 = f0_user(33); 
	fC027 = f0_user(34); 
	fC028 = f0_user(35); 
else
	fmu01 = 1; 
	fmu02 = 0; 
	fmu03 = 1; 
	fmu04 = 0; 
	fmu05 = 0; 
	fmu06 = 0; 
	fmu07 = 0; 
	fC01 = 0; 
	fC02 = 0; 
	fC03 = 0; 
	fC04 = 0; 
	fC05 = 0; 
	fC06 = 0; 
	fC07 = 0; 
	fC08 = 0; 
	fC09 = 0; 
	fC010 = 0; 
	fC011 = 0; 
	fC012 = 0; 
	fC013 = 0; 
	fC014 = 0; 
	fC015 = 0; 
	fC016 = 0; 
	fC017 = 0; 
	fC018 = 0; 
	fC019 = 0; 
	fC020 = 0; 
	fC021 = 0; 
	fC022 = 0; 
	fC023 = 0; 
	fC024 = 0; 
	fC025 = 0; 
	fC026 = 0; 
	fC027 = 0; 
	fC028 = 0; 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = -(C_1_5*Omega^2*kon1 - Omega*koff1*mu_2 + Omega^2*kon1*mu_1*mu_5)/Omega;
xdot(2) = (C_1_5*Omega^2*kon1 - Omega*koff1*mu_2 + Omega^2*kon1*mu_1*mu_5)/Omega;
xdot(3) = -(C_3_6*Omega^2*kon2 - Omega*koff2*mu_4 + Omega^2*kon2*mu_3*mu_6)/Omega;
xdot(4) = (C_3_6*Omega^2*kon2 - Omega*koff2*mu_4 + Omega^2*kon2*mu_3*mu_6)/Omega;
xdot(5) = -(Omega*betam1*mu_5 - 4*beta1 - Omega*koff1*mu_2 + C_1_5*Omega^2*kon1 + Omega^2*kon1*mu_1*mu_5)/Omega;
xdot(6) = -(Omega*betam2*mu_6 - 12*beta2 - Omega*koff2*mu_4 + C_3_6*Omega^2*kon2 + Omega^2*kon2*mu_3*mu_6)/Omega;
xdot(7) = (Omega*mu_2*r1 - Omega*mu_7*myalpha + Omega*mu_4*r2)/Omega;
xdot(8) = (Omega*koff1*mu_2 + 2*C_1_2*Omega^2*koff1 + C_1_5*Omega^2*kon1 - 2*C_1_1*Omega^3*kon1*mu_5 - 2*C_1_5*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_5)/Omega^2;
xdot(9) = -(Omega*koff1*mu_2 - kon1*(C_1_2*Omega^3*mu_5 - Omega*mu_2*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_2_5*Omega^2 + Omega^2*mu_2*mu_5) - Omega*mu_5*(C_1_2*Omega^2 + Omega^2*mu_1*mu_2) + C_1_5*Omega^3*mu_2 + C_2_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_2*mu_5) + C_1_2*Omega^2*koff1 - C_2_2*Omega^2*koff1 + C_1_5*Omega^2*kon1 - C_1_1*Omega^3*kon1*mu_5 - C_1_5*Omega^3*kon1*mu_1 + C_1_2*Omega^3*kon1*mu_5 + C_2_5*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_5)/Omega^2;
xdot(10) = (kon1*(C_1_3*Omega^3*mu_5 - Omega*mu_3*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) - Omega*mu_5*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_5*Omega^3*mu_3 + C_3_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_5) + kon2*(C_1_3*Omega^3*mu_6 - Omega*mu_3*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) - Omega*mu_1*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_6*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_6*Omega^3*mu_3 + C_3_6*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_6) + C_1_4*Omega^2*koff2 + C_2_3*Omega^2*koff1 - C_1_3*Omega^3*kon1*mu_5 - C_1_3*Omega^3*kon2*mu_6 - C_1_6*Omega^3*kon2*mu_3 - C_3_5*Omega^3*kon1*mu_1)/Omega^2;
xdot(11) = -(kon2*(C_1_3*Omega^3*mu_6 - Omega*mu_3*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) - Omega*mu_1*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_6*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_6*Omega^3*mu_3 + C_3_6*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_6) - kon1*(C_1_4*Omega^3*mu_5 - Omega*mu_4*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_4_5*Omega^2 + Omega^2*mu_4*mu_5) - Omega*mu_5*(C_1_4*Omega^2 + Omega^2*mu_1*mu_4) + C_1_5*Omega^3*mu_4 + C_4_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_4*mu_5) + C_1_4*Omega^2*koff2 - C_2_4*Omega^2*koff1 + C_1_4*Omega^3*kon1*mu_5 - C_1_3*Omega^3*kon2*mu_6 - C_1_6*Omega^3*kon2*mu_3 + C_4_5*Omega^3*kon1*mu_1)/Omega^2;
xdot(12) = (Omega*koff1*mu_2 - C_1_5*Omega^2*betam1 + C_1_2*Omega^2*koff1 + C_2_5*Omega^2*koff1 + C_1_5*Omega^2*kon1 - C_1_1*Omega^3*kon1*mu_5 - C_1_5*Omega^3*kon1*mu_1 - C_1_5*Omega^3*kon1*mu_5 - C_5_5*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_5)/Omega^2;
xdot(13) = -(C_1_6*Omega^2*betam2 - kon1*(C_1_5*Omega^3*mu_6 - Omega*mu_5*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) - Omega*mu_1*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) - Omega*mu_6*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_6*Omega^3*mu_5 + C_5_6*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_6) - kon2*(C_1_3*Omega^3*mu_6 - Omega*mu_3*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) - Omega*mu_1*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_6*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_6*Omega^3*mu_3 + C_3_6*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_6) - C_1_4*Omega^2*koff2 - C_2_6*Omega^2*koff1 + C_1_3*Omega^3*kon2*mu_6 + C_1_6*Omega^3*kon2*mu_3 + C_1_6*Omega^3*kon1*mu_5 + C_5_6*Omega^3*kon1*mu_1)/Omega^2;
xdot(14) = (kon1*(C_1_5*Omega^3*mu_7 - Omega*mu_5*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) - Omega*mu_1*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) - Omega*mu_7*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_7*Omega^3*mu_5 + C_5_7*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_7) + C_2_7*Omega^2*koff1 - C_1_7*Omega^2*myalpha + C_1_2*Omega^2*r1 + C_1_4*Omega^2*r2 - C_1_7*Omega^3*kon1*mu_5 - C_5_7*Omega^3*kon1*mu_1)/Omega^2;
xdot(15) = (Omega*koff1*mu_2 - 2*kon1*(C_1_2*Omega^3*mu_5 - Omega*mu_2*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_2_5*Omega^2 + Omega^2*mu_2*mu_5) - Omega*mu_5*(C_1_2*Omega^2 + Omega^2*mu_1*mu_2) + C_1_5*Omega^3*mu_2 + C_2_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_2*mu_5) - 2*C_2_2*Omega^2*koff1 + C_1_5*Omega^2*kon1 + 2*C_1_2*Omega^3*kon1*mu_5 + 2*C_2_5*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_5)/Omega^2;
xdot(16) = -(kon1*(C_1_3*Omega^3*mu_5 - Omega*mu_3*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) - Omega*mu_5*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_5*Omega^3*mu_3 + C_3_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_5) - kon2*(C_2_3*Omega^3*mu_6 - Omega*mu_3*(C_2_6*Omega^2 + Omega^2*mu_2*mu_6) - Omega*mu_2*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_6*(C_2_3*Omega^2 + Omega^2*mu_2*mu_3) + C_2_6*Omega^3*mu_3 + C_3_6*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_3*mu_6) + C_2_3*Omega^2*koff1 - C_2_4*Omega^2*koff2 - C_1_3*Omega^3*kon1*mu_5 + C_2_3*Omega^3*kon2*mu_6 + C_2_6*Omega^3*kon2*mu_3 - C_3_5*Omega^3*kon1*mu_1)/Omega^2;
xdot(17) = -(kon1*(C_1_4*Omega^3*mu_5 - Omega*mu_4*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_4_5*Omega^2 + Omega^2*mu_4*mu_5) - Omega*mu_5*(C_1_4*Omega^2 + Omega^2*mu_1*mu_4) + C_1_5*Omega^3*mu_4 + C_4_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_4*mu_5) + kon2*(C_2_3*Omega^3*mu_6 - Omega*mu_3*(C_2_6*Omega^2 + Omega^2*mu_2*mu_6) - Omega*mu_2*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_6*(C_2_3*Omega^2 + Omega^2*mu_2*mu_3) + C_2_6*Omega^3*mu_3 + C_3_6*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_3*mu_6) + C_2_4*Omega^2*koff1 + C_2_4*Omega^2*koff2 - C_1_4*Omega^3*kon1*mu_5 - C_2_3*Omega^3*kon2*mu_6 - C_2_6*Omega^3*kon2*mu_3 - C_4_5*Omega^3*kon1*mu_1)/Omega^2;
xdot(18) = -(Omega*koff1*mu_2 - kon1*(C_1_2*Omega^3*mu_5 - Omega*mu_2*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_2_5*Omega^2 + Omega^2*mu_2*mu_5) - Omega*mu_5*(C_1_2*Omega^2 + Omega^2*mu_1*mu_2) + C_1_5*Omega^3*mu_2 + C_2_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_2*mu_5) + C_2_5*Omega^2*betam1 - C_2_2*Omega^2*koff1 + C_2_5*Omega^2*koff1 + C_1_5*Omega^2*kon1 + C_1_2*Omega^3*kon1*mu_5 - C_1_5*Omega^3*kon1*mu_5 + C_2_5*Omega^3*kon1*mu_1 - C_5_5*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_5)/Omega^2;
xdot(19) = -(kon1*(C_1_5*Omega^3*mu_6 - Omega*mu_5*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) - Omega*mu_1*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) - Omega*mu_6*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_6*Omega^3*mu_5 + C_5_6*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_6) - kon2*(C_2_3*Omega^3*mu_6 - Omega*mu_3*(C_2_6*Omega^2 + Omega^2*mu_2*mu_6) - Omega*mu_2*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_6*(C_2_3*Omega^2 + Omega^2*mu_2*mu_3) + C_2_6*Omega^3*mu_3 + C_3_6*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_3*mu_6) + C_2_6*Omega^2*betam2 - C_2_4*Omega^2*koff2 + C_2_6*Omega^2*koff1 - C_1_6*Omega^3*kon1*mu_5 + C_2_3*Omega^3*kon2*mu_6 + C_2_6*Omega^3*kon2*mu_3 - C_5_6*Omega^3*kon1*mu_1)/Omega^2;
xdot(20) = (C_2_2*Omega^2*r1 - C_2_7*Omega^2*koff1 - C_2_7*Omega^2*myalpha - kon1*(C_1_5*Omega^3*mu_7 - Omega*mu_5*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) - Omega*mu_1*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) - Omega*mu_7*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_7*Omega^3*mu_5 + C_5_7*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_7) + C_2_4*Omega^2*r2 + C_1_7*Omega^3*kon1*mu_5 + C_5_7*Omega^3*kon1*mu_1)/Omega^2;
xdot(21) = (Omega*koff2*mu_4 + 2*C_3_4*Omega^2*koff2 + C_3_6*Omega^2*kon2 - 2*C_3_3*Omega^3*kon2*mu_6 - 2*C_3_6*Omega^3*kon2*mu_3 + Omega^2*kon2*mu_3*mu_6)/Omega^2;
xdot(22) = -(Omega*koff2*mu_4 - kon2*(C_3_4*Omega^3*mu_6 - Omega*mu_4*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_3*(C_4_6*Omega^2 + Omega^2*mu_4*mu_6) - Omega*mu_6*(C_3_4*Omega^2 + Omega^2*mu_3*mu_4) + C_3_6*Omega^3*mu_4 + C_4_6*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_4*mu_6) + C_3_4*Omega^2*koff2 - C_4_4*Omega^2*koff2 + C_3_6*Omega^2*kon2 - C_3_3*Omega^3*kon2*mu_6 - C_3_6*Omega^3*kon2*mu_3 + C_3_4*Omega^3*kon2*mu_6 + C_4_6*Omega^3*kon2*mu_3 + Omega^2*kon2*mu_3*mu_6)/Omega^2;
xdot(23) = -(C_3_5*Omega^2*betam1 - kon2*(C_3_5*Omega^3*mu_6 - Omega*mu_5*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_3*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) - Omega*mu_6*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_6*Omega^3*mu_5 + C_5_6*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_6) - kon1*(C_1_3*Omega^3*mu_5 - Omega*mu_3*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) - Omega*mu_5*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_5*Omega^3*mu_3 + C_3_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_5) - C_2_3*Omega^2*koff1 - C_4_5*Omega^2*koff2 + C_1_3*Omega^3*kon1*mu_5 + C_3_5*Omega^3*kon1*mu_1 + C_3_5*Omega^3*kon2*mu_6 + C_5_6*Omega^3*kon2*mu_3)/Omega^2;
xdot(24) = (Omega*koff2*mu_4 - C_3_6*Omega^2*betam2 + C_3_4*Omega^2*koff2 + C_4_6*Omega^2*koff2 + C_3_6*Omega^2*kon2 - C_3_3*Omega^3*kon2*mu_6 - C_3_6*Omega^3*kon2*mu_3 - C_3_6*Omega^3*kon2*mu_6 - C_6_6*Omega^3*kon2*mu_3 + Omega^2*kon2*mu_3*mu_6)/Omega^2;
xdot(25) = (kon2*(C_3_6*Omega^3*mu_7 - Omega*mu_6*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) - Omega*mu_3*(C_6_7*Omega^2 + Omega^2*mu_6*mu_7) - Omega*mu_7*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) + C_3_7*Omega^3*mu_6 + C_6_7*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_6*mu_7) + C_4_7*Omega^2*koff2 - C_3_7*Omega^2*myalpha + C_2_3*Omega^2*r1 + C_3_4*Omega^2*r2 - C_3_7*Omega^3*kon2*mu_6 - C_6_7*Omega^3*kon2*mu_3)/Omega^2;
xdot(26) = (Omega*koff2*mu_4 - 2*kon2*(C_3_4*Omega^3*mu_6 - Omega*mu_4*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_3*(C_4_6*Omega^2 + Omega^2*mu_4*mu_6) - Omega*mu_6*(C_3_4*Omega^2 + Omega^2*mu_3*mu_4) + C_3_6*Omega^3*mu_4 + C_4_6*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_4*mu_6) - 2*C_4_4*Omega^2*koff2 + C_3_6*Omega^2*kon2 + 2*C_3_4*Omega^3*kon2*mu_6 + 2*C_4_6*Omega^3*kon2*mu_3 + Omega^2*kon2*mu_3*mu_6)/Omega^2;
xdot(27) = -(kon2*(C_3_5*Omega^3*mu_6 - Omega*mu_5*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_3*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) - Omega*mu_6*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_6*Omega^3*mu_5 + C_5_6*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_6) - kon1*(C_1_4*Omega^3*mu_5 - Omega*mu_4*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) - Omega*mu_1*(C_4_5*Omega^2 + Omega^2*mu_4*mu_5) - Omega*mu_5*(C_1_4*Omega^2 + Omega^2*mu_1*mu_4) + C_1_5*Omega^3*mu_4 + C_4_5*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_4*mu_5) + C_4_5*Omega^2*betam1 - C_2_4*Omega^2*koff1 + C_4_5*Omega^2*koff2 + C_1_4*Omega^3*kon1*mu_5 - C_3_5*Omega^3*kon2*mu_6 + C_4_5*Omega^3*kon1*mu_1 - C_5_6*Omega^3*kon2*mu_3)/Omega^2;
xdot(28) = -(Omega*koff2*mu_4 - kon2*(C_3_4*Omega^3*mu_6 - Omega*mu_4*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_3*(C_4_6*Omega^2 + Omega^2*mu_4*mu_6) - Omega*mu_6*(C_3_4*Omega^2 + Omega^2*mu_3*mu_4) + C_3_6*Omega^3*mu_4 + C_4_6*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_4*mu_6) + C_4_6*Omega^2*betam2 - C_4_4*Omega^2*koff2 + C_4_6*Omega^2*koff2 + C_3_6*Omega^2*kon2 + C_3_4*Omega^3*kon2*mu_6 - C_3_6*Omega^3*kon2*mu_6 + C_4_6*Omega^3*kon2*mu_3 - C_6_6*Omega^3*kon2*mu_3 + Omega^2*kon2*mu_3*mu_6)/Omega^2;
xdot(29) = (C_2_4*Omega^2*r1 - C_4_7*Omega^2*koff2 - C_4_7*Omega^2*myalpha - kon2*(C_3_6*Omega^3*mu_7 - Omega*mu_6*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) - Omega*mu_3*(C_6_7*Omega^2 + Omega^2*mu_6*mu_7) - Omega*mu_7*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) + C_3_7*Omega^3*mu_6 + C_6_7*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_6*mu_7) + C_4_4*Omega^2*r2 + C_3_7*Omega^3*kon2*mu_6 + C_6_7*Omega^3*kon2*mu_3)/Omega^2;
xdot(30) = (16*beta1 + Omega*betam1*mu_5 + Omega*koff1*mu_2 - 2*C_5_5*Omega^2*betam1 + 2*C_2_5*Omega^2*koff1 + C_1_5*Omega^2*kon1 - 2*C_1_5*Omega^3*kon1*mu_5 - 2*C_5_5*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_5)/Omega^2;
xdot(31) = -(C_5_6*Omega^2*betam1 - kon2*(C_3_5*Omega^3*mu_6 - Omega*mu_5*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) - Omega*mu_3*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) - Omega*mu_6*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_6*Omega^3*mu_5 + C_5_6*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_6) - kon1*(C_1_5*Omega^3*mu_6 - Omega*mu_5*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) - Omega*mu_1*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) - Omega*mu_6*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_6*Omega^3*mu_5 + C_5_6*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_6) + C_5_6*Omega^2*betam2 - C_2_6*Omega^2*koff1 - C_4_5*Omega^2*koff2 + C_1_6*Omega^3*kon1*mu_5 + C_3_5*Omega^3*kon2*mu_6 + C_5_6*Omega^3*kon1*mu_1 + C_5_6*Omega^3*kon2*mu_3)/Omega^2;
xdot(32) = (kon1*(C_1_5*Omega^3*mu_7 - Omega*mu_5*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) - Omega*mu_1*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) - Omega*mu_7*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_7*Omega^3*mu_5 + C_5_7*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_7) - C_5_7*Omega^2*betam1 + C_2_7*Omega^2*koff1 - C_5_7*Omega^2*myalpha + C_2_5*Omega^2*r1 + C_4_5*Omega^2*r2 - C_1_7*Omega^3*kon1*mu_5 - C_5_7*Omega^3*kon1*mu_1)/Omega^2;
xdot(33) = (144*beta2 + Omega*betam2*mu_6 + Omega*koff2*mu_4 - 2*C_6_6*Omega^2*betam2 + 2*C_4_6*Omega^2*koff2 + C_3_6*Omega^2*kon2 - 2*C_3_6*Omega^3*kon2*mu_6 - 2*C_6_6*Omega^3*kon2*mu_3 + Omega^2*kon2*mu_3*mu_6)/Omega^2;
xdot(34) = (kon2*(C_3_6*Omega^3*mu_7 - Omega*mu_6*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) - Omega*mu_3*(C_6_7*Omega^2 + Omega^2*mu_6*mu_7) - Omega*mu_7*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) + C_3_7*Omega^3*mu_6 + C_6_7*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_6*mu_7) - C_6_7*Omega^2*betam2 + C_4_7*Omega^2*koff2 - C_6_7*Omega^2*myalpha + C_2_6*Omega^2*r1 + C_4_6*Omega^2*r2 - C_3_7*Omega^3*kon2*mu_6 - C_6_7*Omega^3*kon2*mu_3)/Omega^2;
xdot(35) = (Omega*mu_7*myalpha + Omega*mu_2*r1 + Omega*mu_4*r2 - 2*C_7_7*Omega^2*myalpha + 2*C_2_7*Omega^2*r1 + 2*C_4_7*Omega^2*r2)/Omega^2;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01 - fmu01*(indmu1 - 1))/Omega;
x0(2) = (indmu2*kmu02 - fmu02*(indmu2 - 1))/Omega;
x0(3) = (indmu3*kmu03 - fmu03*(indmu3 - 1))/Omega;
x0(4) = (indmu4*kmu04 - fmu04*(indmu4 - 1))/Omega;
x0(5) = (indmu5*kmu05 - fmu05*(indmu5 - 1))/Omega;
x0(6) = (indmu6*kmu06 - fmu06*(indmu6 - 1))/Omega;
x0(7) = (indmu7*kmu07 - fmu07*(indmu7 - 1))/Omega;
x0(8) = (indC1*kC01 - fC01*(indC1 - 1))/Omega^2;
x0(9) = (indC2*kC02 - fC02*(indC2 - 1))/Omega^2;
x0(10) = (indC3*kC03 - fC03*(indC3 - 1))/Omega^2;
x0(11) = (indC4*kC04 - fC04*(indC4 - 1))/Omega^2;
x0(12) = (indC5*kC05 - fC05*(indC5 - 1))/Omega^2;
x0(13) = (indC6*kC06 - fC06*(indC6 - 1))/Omega^2;
x0(14) = (indC7*kC07 - fC07*(indC7 - 1))/Omega^2;
x0(15) = (indC8*kC08 - fC08*(indC8 - 1))/Omega^2;
x0(16) = (indC9*kC09 - fC09*(indC9 - 1))/Omega^2;
x0(17) = (indC10*kC010 - fC010*(indC10 - 1))/Omega^2;
x0(18) = (indC11*kC011 - fC011*(indC11 - 1))/Omega^2;
x0(19) = (indC12*kC012 - fC012*(indC12 - 1))/Omega^2;
x0(20) = (indC13*kC013 - fC013*(indC13 - 1))/Omega^2;
x0(21) = (indC14*kC014 - fC014*(indC14 - 1))/Omega^2;
x0(22) = (indC15*kC015 - fC015*(indC15 - 1))/Omega^2;
x0(23) = (indC16*kC016 - fC016*(indC16 - 1))/Omega^2;
x0(24) = (indC17*kC017 - fC017*(indC17 - 1))/Omega^2;
x0(25) = (indC18*kC018 - fC018*(indC18 - 1))/Omega^2;
x0(26) = (indC19*kC019 - fC019*(indC19 - 1))/Omega^2;
x0(27) = (indC20*kC020 - fC020*(indC20 - 1))/Omega^2;
x0(28) = (indC21*kC021 - fC021*(indC21 - 1))/Omega^2;
x0(29) = (indC22*kC022 - fC022*(indC22 - 1))/Omega^2;
x0(30) = (indC23*kC023 - fC023*(indC23 - 1))/Omega^2;
x0(31) = (indC24*kC024 - fC024*(indC24 - 1))/Omega^2;
x0(32) = (indC25*kC025 - fC025*(indC25 - 1))/Omega^2;
x0(33) = (indC26*kC026 - fC026*(indC26 - 1))/Omega^2;
x0(34) = (indC27*kC027 - fC027*(indC27 - 1))/Omega^2;
x0(35) = (indC28*kC028 - fC028*(indC28 - 1))/Omega^2;

% OBSERVABLES

y = sym(zeros(2,1));

y(1) = mu_7;
y(2) = C_7_7;

% SYSTEM STRUCT

model.sym.nmx = 0;
model.sym.x = x;
model.sym.u = u;
model.sym.xdot = xdot;
model.sym.p = p;
model.sym.k = k;
model.sym.x0 = x0;
model.sym.y = y;
% Additional fields for the prespecified length of kappa
model.sym.nk1 = 1;
end