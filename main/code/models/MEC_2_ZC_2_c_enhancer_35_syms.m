% function [model] = MEC_2_ZC_2_c_enhancer_35_syms(f0_user)
function [model] = MEC_2_ZC_2_c_enhancer_35_syms(varargin)

% CVODES OPTIONS

model.atol = 1e-8;
model.rtol = 1e-8;
model.maxsteps = 1e4;

% STATES

syms mu_1 mu_2 mu_3 mu_4 mu_5 mu_6 mu_7 mu_8 mu_9 mu_10 mu_11 C_1_1 C_1_2 C_1_3 C_1_4 C_1_5 C_1_6 C_1_7 C_1_8 C_1_9 C_1_10 C_1_11 C_2_2 C_2_3 C_2_4 C_2_5 C_2_6 C_2_7 C_2_8 C_2_9 C_2_10 C_2_11 C_3_3 C_3_4 C_3_5 C_3_6 C_3_7 C_3_8 C_3_9 C_3_10 C_3_11 C_4_4 C_4_5 C_4_6 C_4_7 C_4_8 C_4_9 C_4_10 C_4_11 C_5_5 C_5_6 C_5_7 C_5_8 C_5_9 C_5_10 C_5_11 C_6_6 C_6_7 C_6_8 C_6_9 C_6_10 C_6_11 C_7_7 C_7_8 C_7_9 C_7_10 C_7_11 C_8_8 C_8_9 C_8_10 C_8_11 C_9_9 C_9_10 C_9_11 C_10_10 C_10_11 C_11_11

x = [
mu_1, mu_2, mu_3, mu_4, mu_5, mu_6, mu_7, mu_8, mu_9, mu_10, mu_11, C_1_1, C_1_2, C_1_3, C_1_4, C_1_5, C_1_6, C_1_7, C_1_8, C_1_9, C_1_10, C_1_11, C_2_2, C_2_3, C_2_4, C_2_5, C_2_6, C_2_7, C_2_8, C_2_9, C_2_10, C_2_11, C_3_3, C_3_4, C_3_5, C_3_6, C_3_7, C_3_8, C_3_9, C_3_10, C_3_11, C_4_4, C_4_5, C_4_6, C_4_7, C_4_8, C_4_9, C_4_10, C_4_11, C_5_5, C_5_6, C_5_7, C_5_8, C_5_9, C_5_10, C_5_11, C_6_6, C_6_7, C_6_8, C_6_9, C_6_10, C_6_11, C_7_7, C_7_8, C_7_9, C_7_10, C_7_11, C_8_8, C_8_9, C_8_10, C_8_11, C_9_9, C_9_10, C_9_11, C_10_10, C_10_11, C_11_11 ...
];

% PARAMETERS

syms kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 

% KAPPA (constant parameters)

syms Omega indmu1 indmu2 indmu3 indmu4 indmu5 indmu6 indmu7 indmu8 indmu9 indmu10 indmu11 indC1 indC2 indC3 indC4 indC5 indC6 indC7 indC8 indC9 indC10 indC11 indC12 indC13 indC14 indC15 indC16 indC17 indC18 indC19 indC20 indC21 indC22 indC23 indC24 indC25 indC26 indC27 indC28 indC29 indC30 indC31 indC32 indC33 indC34 indC35 indC36 indC37 indC38 indC39 indC40 indC41 indC42 indC43 indC44 indC45 indC46 indC47 indC48 indC49 indC50 indC51 indC52 indC53 indC54 indC55 indC56 indC57 indC58 indC59 indC60 indC61 indC62 indC63 indC64 indC65 indC66 kmu01 kmu02 kmu03 kmu04 kmu05 kmu06 kmu07 kmu08 kmu09 kmu010 kmu011 kC01 kC02 kC03 kC04 kC05 kC06 kC07 kC08 kC09 kC010 kC011 kC012 kC013 kC014 kC015 kC016 kC017 kC018 kC019 kC020 kC021 kC022 kC023 kC024 kC025 kC026 kC027 kC028 kC029 kC030 kC031 kC032 kC033 kC034 kC035 kC036 kC037 kC038 kC039 kC040 kC041 kC042 kC043 kC044 kC045 kC046 kC047 kC048 kC049 kC050 kC051 kC052 kC053 kC054 kC055 kC056 kC057 kC058 kC059 kC060 kC061 kC062 kC063 kC064 kC065 kC066 

syms t

p = [kon1,koff1,r1,kon2,koff2,r2,myalpha,beta1,betam1,beta2,betam2];

k = [Omega,indmu1,indmu2,indmu3,indmu4,indmu5,indmu6,indmu7,indmu8,indmu9,indmu10,indmu11,indC1,indC2,indC3,indC4,indC5,indC6,indC7,indC8,indC9,indC10,indC11,indC12,indC13,indC14,indC15,indC16,indC17,indC18,indC19,indC20,indC21,indC22,indC23,indC24,indC25,indC26,indC27,indC28,indC29,indC30,indC31,indC32,indC33,indC34,indC35,indC36,indC37,indC38,indC39,indC40,indC41,indC42,indC43,indC44,indC45,indC46,indC47,indC48,indC49,indC50,indC51,indC52,indC53,indC54,indC55,indC56,indC57,indC58,indC59,indC60,indC61,indC62,indC63,indC64,indC65,indC66,kmu01,kmu02,kmu03,kmu04,kmu05,kmu06,kmu07,kmu08,kmu09,kmu010,kmu011,kC01,kC02,kC03,kC04,kC05,kC06,kC07,kC08,kC09,kC010,kC011,kC012,kC013,kC014,kC015,kC016,kC017,kC018,kC019,kC020,kC021,kC022,kC023,kC024,kC025,kC026,kC027,kC028,kC029,kC030,kC031,kC032,kC033,kC034,kC035,kC036,kC037,kC038,kC039,kC040,kC041,kC042,kC043,kC044,kC045,kC046,kC047,kC048,kC049,kC050,kC051,kC052,kC053,kC054,kC055,kC056,kC057,kC058,kC059,kC060,kC061,kC062,kC063,kC064,kC065,kC066];

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
	fmu08 = f0_user(8); 
	fmu09 = f0_user(9); 
	fmu010 = f0_user(10); 
	fmu011 = f0_user(11); 
	fC01 = f0_user(12); 
	fC02 = f0_user(13); 
	fC03 = f0_user(14); 
	fC04 = f0_user(15); 
	fC05 = f0_user(16); 
	fC06 = f0_user(17); 
	fC07 = f0_user(18); 
	fC08 = f0_user(19); 
	fC09 = f0_user(20); 
	fC010 = f0_user(21); 
	fC011 = f0_user(22); 
	fC012 = f0_user(23); 
	fC013 = f0_user(24); 
	fC014 = f0_user(25); 
	fC015 = f0_user(26); 
	fC016 = f0_user(27); 
	fC017 = f0_user(28); 
	fC018 = f0_user(29); 
	fC019 = f0_user(30); 
	fC020 = f0_user(31); 
	fC021 = f0_user(32); 
	fC022 = f0_user(33); 
	fC023 = f0_user(34); 
	fC024 = f0_user(35); 
	fC025 = f0_user(36); 
	fC026 = f0_user(37); 
	fC027 = f0_user(38); 
	fC028 = f0_user(39); 
	fC029 = f0_user(40); 
	fC030 = f0_user(41); 
	fC031 = f0_user(42); 
	fC032 = f0_user(43); 
	fC033 = f0_user(44); 
	fC034 = f0_user(45); 
	fC035 = f0_user(46); 
	fC036 = f0_user(47); 
	fC037 = f0_user(48); 
	fC038 = f0_user(49); 
	fC039 = f0_user(50); 
	fC040 = f0_user(51); 
	fC041 = f0_user(52); 
	fC042 = f0_user(53); 
	fC043 = f0_user(54); 
	fC044 = f0_user(55); 
	fC045 = f0_user(56); 
	fC046 = f0_user(57); 
	fC047 = f0_user(58); 
	fC048 = f0_user(59); 
	fC049 = f0_user(60); 
	fC050 = f0_user(61); 
	fC051 = f0_user(62); 
	fC052 = f0_user(63); 
	fC053 = f0_user(64); 
	fC054 = f0_user(65); 
	fC055 = f0_user(66); 
	fC056 = f0_user(67); 
	fC057 = f0_user(68); 
	fC058 = f0_user(69); 
	fC059 = f0_user(70); 
	fC060 = f0_user(71); 
	fC061 = f0_user(72); 
	fC062 = f0_user(73); 
	fC063 = f0_user(74); 
	fC064 = f0_user(75); 
	fC065 = f0_user(76); 
	fC066 = f0_user(77); 
else
	fmu01 = 1; 
	fmu02 = 0; 
	fmu03 = 1; 
	fmu04 = 0; 
	fmu05 = 1; 
	fmu06 = 0; 
	fmu07 = 1; 
	fmu08 = 0; 
	fmu09 = 0; 
	fmu010 = 0; 
	fmu011 = 0; 
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
	fC029 = 0; 
	fC030 = 0; 
	fC031 = 0; 
	fC032 = 0; 
	fC033 = 0; 
	fC034 = 0; 
	fC035 = 0; 
	fC036 = 0; 
	fC037 = 0; 
	fC038 = 0; 
	fC039 = 0; 
	fC040 = 0; 
	fC041 = 0; 
	fC042 = 0; 
	fC043 = 0; 
	fC044 = 0; 
	fC045 = 0; 
	fC046 = 0; 
	fC047 = 0; 
	fC048 = 0; 
	fC049 = 0; 
	fC050 = 0; 
	fC051 = 0; 
	fC052 = 0; 
	fC053 = 0; 
	fC054 = 0; 
	fC055 = 0; 
	fC056 = 0; 
	fC057 = 0; 
	fC058 = 0; 
	fC059 = 0; 
	fC060 = 0; 
	fC061 = 0; 
	fC062 = 0; 
	fC063 = 0; 
	fC064 = 0; 
	fC065 = 0; 
	fC066 = 0; 
end
% INPUT 

u = sym.empty(0,0);

% SYSTEM EQUATIONS

xdot = sym(zeros(size(x)));

xdot(1) = -(C_1_9*Omega^2*kon1 - Omega*koff1*mu_2 + Omega^2*kon1*mu_1*mu_9)/Omega;
xdot(2) = (C_1_9*Omega^2*kon1 - Omega*koff1*mu_2 + Omega^2*kon1*mu_1*mu_9)/Omega;
xdot(3) = -(C_3_9*Omega^2*kon1 - Omega*koff1*mu_4 + Omega^2*kon1*mu_3*mu_9)/Omega;
xdot(4) = (C_3_9*Omega^2*kon1 - Omega*koff1*mu_4 + Omega^2*kon1*mu_3*mu_9)/Omega;
xdot(5) = -(C_5_10*Omega^2*kon2 - Omega*koff2*mu_6 + Omega^2*kon2*mu_5*mu_10)/Omega;
xdot(6) = (C_5_10*Omega^2*kon2 - Omega*koff2*mu_6 + Omega^2*kon2*mu_5*mu_10)/Omega;
xdot(7) = -(C_7_10*Omega^2*kon2 - Omega*koff2*mu_8 + Omega^2*kon2*mu_7*mu_10)/Omega;
xdot(8) = (C_7_10*Omega^2*kon2 - Omega*koff2*mu_8 + Omega^2*kon2*mu_7*mu_10)/Omega;
xdot(9) = -(Omega*betam1*mu_9 - 4*beta1 - Omega*koff1*mu_2 - Omega*koff1*mu_4 + C_1_9*Omega^2*kon1 + C_3_9*Omega^2*kon1 + Omega^2*kon1*mu_1*mu_9 + Omega^2*kon1*mu_3*mu_9)/Omega;
xdot(10) = -(Omega*betam2*mu_10 - 12*beta2 - Omega*koff2*mu_6 - Omega*koff2*mu_8 + C_5_10*Omega^2*kon2 + C_7_10*Omega^2*kon2 + Omega^2*kon2*mu_5*mu_10 + Omega^2*kon2*mu_7*mu_10)/Omega;
xdot(11) = (Omega*mu_2*r1 - Omega*mu_11*myalpha + Omega*mu_4*r1 + Omega*mu_6*r2 + Omega*mu_8*r2)/Omega;
xdot(12) = (Omega*koff1*mu_2 + 2*C_1_2*Omega^2*koff1 + C_1_9*Omega^2*kon1 - 2*C_1_1*Omega^3*kon1*mu_9 - 2*C_1_9*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_9)/Omega^2;
xdot(13) = -(Omega*koff1*mu_2 - kon1*(C_1_2*Omega^3*mu_9 - Omega*mu_2*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_2_9*Omega^2 + Omega^2*mu_2*mu_9) - Omega*mu_9*(C_1_2*Omega^2 + Omega^2*mu_1*mu_2) + C_1_9*Omega^3*mu_2 + C_2_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_2*mu_9) + C_1_2*Omega^2*koff1 - C_2_2*Omega^2*koff1 + C_1_9*Omega^2*kon1 - C_1_1*Omega^3*kon1*mu_9 - C_1_9*Omega^3*kon1*mu_1 + C_1_2*Omega^3*kon1*mu_9 + C_2_9*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_9)/Omega^2;
xdot(14) = (2*kon1*(C_1_3*Omega^3*mu_9 - Omega*mu_3*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_9) + C_1_4*Omega^2*koff1 + C_2_3*Omega^2*koff1 - 2*C_1_3*Omega^3*kon1*mu_9 - C_1_9*Omega^3*kon1*mu_3 - C_3_9*Omega^3*kon1*mu_1)/Omega^2;
xdot(15) = -(kon1*(C_1_3*Omega^3*mu_9 - Omega*mu_3*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_9) - kon1*(C_1_4*Omega^3*mu_9 - Omega*mu_4*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_4_9*Omega^2 + Omega^2*mu_4*mu_9) - Omega*mu_9*(C_1_4*Omega^2 + Omega^2*mu_1*mu_4) + C_1_9*Omega^3*mu_4 + C_4_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_4*mu_9) + C_1_4*Omega^2*koff1 - C_2_4*Omega^2*koff1 - C_1_3*Omega^3*kon1*mu_9 - C_1_9*Omega^3*kon1*mu_3 + C_1_4*Omega^3*kon1*mu_9 + C_4_9*Omega^3*kon1*mu_1)/Omega^2;
xdot(16) = (kon1*(C_1_5*Omega^3*mu_9 - Omega*mu_5*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) - Omega*mu_9*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_9*Omega^3*mu_5 + C_5_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_9) + kon2*(C_1_5*Omega^3*mu_10 - Omega*mu_5*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_10) + C_1_6*Omega^2*koff2 + C_2_5*Omega^2*koff1 - C_1_5*Omega^3*kon1*mu_9 - C_1_5*Omega^3*kon2*mu_10 - C_5_9*Omega^3*kon1*mu_1 - C_1_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(17) = (kon1*(C_1_6*Omega^3*mu_9 - Omega*mu_6*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_6_9*Omega^2 + Omega^2*mu_6*mu_9) - Omega*mu_9*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) + C_1_9*Omega^3*mu_6 + C_6_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_6*mu_9) - kon2*(C_1_5*Omega^3*mu_10 - Omega*mu_5*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_10) - C_1_6*Omega^2*koff2 + C_2_6*Omega^2*koff1 - C_1_6*Omega^3*kon1*mu_9 + C_1_5*Omega^3*kon2*mu_10 - C_6_9*Omega^3*kon1*mu_1 + C_1_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(18) = (kon1*(C_1_7*Omega^3*mu_9 - Omega*mu_7*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) - Omega*mu_9*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) + C_1_9*Omega^3*mu_7 + C_7_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_7*mu_9) + kon2*(C_1_7*Omega^3*mu_10 - Omega*mu_7*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) + C_1_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_7*mu_10) + C_1_8*Omega^2*koff2 + C_2_7*Omega^2*koff1 - C_1_7*Omega^3*kon1*mu_9 - C_1_7*Omega^3*kon2*mu_10 - C_7_9*Omega^3*kon1*mu_1 - C_1_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(19) = (kon1*(C_1_8*Omega^3*mu_9 - Omega*mu_8*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_8_9*Omega^2 + Omega^2*mu_8*mu_9) - Omega*mu_9*(C_1_8*Omega^2 + Omega^2*mu_1*mu_8) + C_1_9*Omega^3*mu_8 + C_8_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_8*mu_9) - kon2*(C_1_7*Omega^3*mu_10 - Omega*mu_7*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) + C_1_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_7*mu_10) - C_1_8*Omega^2*koff2 + C_2_8*Omega^2*koff1 - C_1_8*Omega^3*kon1*mu_9 + C_1_7*Omega^3*kon2*mu_10 - C_8_9*Omega^3*kon1*mu_1 + C_1_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(20) = (kon1*(C_1_3*Omega^3*mu_9 - Omega*mu_3*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_9) + Omega*koff1*mu_2 - C_1_9*Omega^2*betam1 + C_1_2*Omega^2*koff1 + C_1_4*Omega^2*koff1 + C_2_9*Omega^2*koff1 + C_1_9*Omega^2*kon1 - C_1_1*Omega^3*kon1*mu_9 - C_1_9*Omega^3*kon1*mu_1 - C_1_3*Omega^3*kon1*mu_9 - C_1_9*Omega^3*kon1*mu_3 - C_1_9*Omega^3*kon1*mu_9 - C_9_9*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_9)/Omega^2;
xdot(21) = -(C_1_10*Omega^2*betam2 - kon2*(C_1_7*Omega^3*mu_10 - Omega*mu_7*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) + C_1_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_7*mu_10) - kon1*(C_1_9*Omega^3*mu_10 - Omega*mu_9*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) + C_1_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_9*mu_10) - kon2*(C_1_5*Omega^3*mu_10 - Omega*mu_5*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_10) - C_1_6*Omega^2*koff2 - C_1_8*Omega^2*koff2 - C_2_10*Omega^2*koff1 + C_1_5*Omega^3*kon2*mu_10 + C_1_7*Omega^3*kon2*mu_10 + C_1_10*Omega^3*kon2*mu_5 + C_1_10*Omega^3*kon2*mu_7 + C_1_10*Omega^3*kon1*mu_9 + C_9_10*Omega^3*kon1*mu_1)/Omega^2;
xdot(22) = (kon1*(C_1_9*Omega^3*mu_11 - Omega*mu_9*(C_1_11*Omega^2 + Omega^2*mu_1*mu_11) - Omega*mu_1*(C_9_11*Omega^2 + Omega^2*mu_9*mu_11) - Omega*mu_11*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) + C_1_11*Omega^3*mu_9 + C_9_11*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_9*mu_11) + C_2_11*Omega^2*koff1 - C_1_11*Omega^2*myalpha + C_1_2*Omega^2*r1 + C_1_4*Omega^2*r1 + C_1_6*Omega^2*r2 + C_1_8*Omega^2*r2 - C_1_11*Omega^3*kon1*mu_9 - C_9_11*Omega^3*kon1*mu_1)/Omega^2;
xdot(23) = (Omega*koff1*mu_2 - 2*kon1*(C_1_2*Omega^3*mu_9 - Omega*mu_2*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_2_9*Omega^2 + Omega^2*mu_2*mu_9) - Omega*mu_9*(C_1_2*Omega^2 + Omega^2*mu_1*mu_2) + C_1_9*Omega^3*mu_2 + C_2_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_2*mu_9) - 2*C_2_2*Omega^2*koff1 + C_1_9*Omega^2*kon1 + 2*C_1_2*Omega^3*kon1*mu_9 + 2*C_2_9*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_9)/Omega^2;
xdot(24) = -(kon1*(C_1_3*Omega^3*mu_9 - Omega*mu_3*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_9) - kon1*(C_2_3*Omega^3*mu_9 - Omega*mu_3*(C_2_9*Omega^2 + Omega^2*mu_2*mu_9) - Omega*mu_2*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_2_3*Omega^2 + Omega^2*mu_2*mu_3) + C_2_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_3*mu_9) + C_2_3*Omega^2*koff1 - C_2_4*Omega^2*koff1 - C_1_3*Omega^3*kon1*mu_9 + C_2_3*Omega^3*kon1*mu_9 + C_2_9*Omega^3*kon1*mu_3 - C_3_9*Omega^3*kon1*mu_1)/Omega^2;
xdot(25) = (C_1_4*Omega^3*kon1*mu_9 - kon1*(C_2_3*Omega^3*mu_9 - Omega*mu_3*(C_2_9*Omega^2 + Omega^2*mu_2*mu_9) - Omega*mu_2*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_2_3*Omega^2 + Omega^2*mu_2*mu_3) + C_2_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_3*mu_9) - 2*C_2_4*Omega^2*koff1 - kon1*(C_1_4*Omega^3*mu_9 - Omega*mu_4*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_4_9*Omega^2 + Omega^2*mu_4*mu_9) - Omega*mu_9*(C_1_4*Omega^2 + Omega^2*mu_1*mu_4) + C_1_9*Omega^3*mu_4 + C_4_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_4*mu_9) + C_2_3*Omega^3*kon1*mu_9 + C_2_9*Omega^3*kon1*mu_3 + C_4_9*Omega^3*kon1*mu_1)/Omega^2;
xdot(26) = -(kon1*(C_1_5*Omega^3*mu_9 - Omega*mu_5*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) - Omega*mu_9*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_9*Omega^3*mu_5 + C_5_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_9) - kon2*(C_2_5*Omega^3*mu_10 - Omega*mu_5*(C_2_10*Omega^2 + Omega^2*mu_2*mu_10) - Omega*mu_2*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_2_5*Omega^2 + Omega^2*mu_2*mu_5) + C_2_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_5*mu_10) + C_2_5*Omega^2*koff1 - C_2_6*Omega^2*koff2 - C_1_5*Omega^3*kon1*mu_9 + C_2_5*Omega^3*kon2*mu_10 - C_5_9*Omega^3*kon1*mu_1 + C_2_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(27) = -(kon1*(C_1_6*Omega^3*mu_9 - Omega*mu_6*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_6_9*Omega^2 + Omega^2*mu_6*mu_9) - Omega*mu_9*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) + C_1_9*Omega^3*mu_6 + C_6_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_6*mu_9) + kon2*(C_2_5*Omega^3*mu_10 - Omega*mu_5*(C_2_10*Omega^2 + Omega^2*mu_2*mu_10) - Omega*mu_2*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_2_5*Omega^2 + Omega^2*mu_2*mu_5) + C_2_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_5*mu_10) + C_2_6*Omega^2*koff1 + C_2_6*Omega^2*koff2 - C_1_6*Omega^3*kon1*mu_9 - C_2_5*Omega^3*kon2*mu_10 - C_6_9*Omega^3*kon1*mu_1 - C_2_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(28) = -(kon1*(C_1_7*Omega^3*mu_9 - Omega*mu_7*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) - Omega*mu_9*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) + C_1_9*Omega^3*mu_7 + C_7_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_7*mu_9) - kon2*(C_2_7*Omega^3*mu_10 - Omega*mu_7*(C_2_10*Omega^2 + Omega^2*mu_2*mu_10) - Omega*mu_2*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_2_7*Omega^2 + Omega^2*mu_2*mu_7) + C_2_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_7*mu_10) + C_2_7*Omega^2*koff1 - C_2_8*Omega^2*koff2 - C_1_7*Omega^3*kon1*mu_9 + C_2_7*Omega^3*kon2*mu_10 - C_7_9*Omega^3*kon1*mu_1 + C_2_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(29) = -(kon1*(C_1_8*Omega^3*mu_9 - Omega*mu_8*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_8_9*Omega^2 + Omega^2*mu_8*mu_9) - Omega*mu_9*(C_1_8*Omega^2 + Omega^2*mu_1*mu_8) + C_1_9*Omega^3*mu_8 + C_8_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_8*mu_9) + kon2*(C_2_7*Omega^3*mu_10 - Omega*mu_7*(C_2_10*Omega^2 + Omega^2*mu_2*mu_10) - Omega*mu_2*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_2_7*Omega^2 + Omega^2*mu_2*mu_7) + C_2_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_7*mu_10) + C_2_8*Omega^2*koff1 + C_2_8*Omega^2*koff2 - C_1_8*Omega^3*kon1*mu_9 - C_2_7*Omega^3*kon2*mu_10 - C_8_9*Omega^3*kon1*mu_1 - C_2_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(30) = -(Omega*koff1*mu_2 - kon1*(C_2_3*Omega^3*mu_9 - Omega*mu_3*(C_2_9*Omega^2 + Omega^2*mu_2*mu_9) - Omega*mu_2*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_2_3*Omega^2 + Omega^2*mu_2*mu_3) + C_2_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_3*mu_9) - kon1*(C_1_2*Omega^3*mu_9 - Omega*mu_2*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_2_9*Omega^2 + Omega^2*mu_2*mu_9) - Omega*mu_9*(C_1_2*Omega^2 + Omega^2*mu_1*mu_2) + C_1_9*Omega^3*mu_2 + C_2_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_2*mu_9) + C_2_9*Omega^2*betam1 - C_2_2*Omega^2*koff1 - C_2_4*Omega^2*koff1 + C_2_9*Omega^2*koff1 + C_1_9*Omega^2*kon1 + C_1_2*Omega^3*kon1*mu_9 - C_1_9*Omega^3*kon1*mu_9 + C_2_9*Omega^3*kon1*mu_1 + C_2_3*Omega^3*kon1*mu_9 + C_2_9*Omega^3*kon1*mu_3 - C_9_9*Omega^3*kon1*mu_1 + Omega^2*kon1*mu_1*mu_9)/Omega^2;
xdot(31) = -(kon1*(C_1_9*Omega^3*mu_10 - Omega*mu_9*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) + C_1_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_9*mu_10) - kon2*(C_2_7*Omega^3*mu_10 - Omega*mu_7*(C_2_10*Omega^2 + Omega^2*mu_2*mu_10) - Omega*mu_2*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_2_7*Omega^2 + Omega^2*mu_2*mu_7) + C_2_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_7*mu_10) - kon2*(C_2_5*Omega^3*mu_10 - Omega*mu_5*(C_2_10*Omega^2 + Omega^2*mu_2*mu_10) - Omega*mu_2*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_2_5*Omega^2 + Omega^2*mu_2*mu_5) + C_2_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_2 + 3*Omega^3*mu_2*mu_5*mu_10) + C_2_10*Omega^2*betam2 - C_2_6*Omega^2*koff2 - C_2_8*Omega^2*koff2 + C_2_10*Omega^2*koff1 + C_2_5*Omega^3*kon2*mu_10 + C_2_7*Omega^3*kon2*mu_10 - C_1_10*Omega^3*kon1*mu_9 + C_2_10*Omega^3*kon2*mu_5 + C_2_10*Omega^3*kon2*mu_7 - C_9_10*Omega^3*kon1*mu_1)/Omega^2;
xdot(32) = (C_2_2*Omega^2*r1 - C_2_11*Omega^2*koff1 - C_2_11*Omega^2*myalpha - kon1*(C_1_9*Omega^3*mu_11 - Omega*mu_9*(C_1_11*Omega^2 + Omega^2*mu_1*mu_11) - Omega*mu_1*(C_9_11*Omega^2 + Omega^2*mu_9*mu_11) - Omega*mu_11*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) + C_1_11*Omega^3*mu_9 + C_9_11*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_9*mu_11) + C_2_4*Omega^2*r1 + C_2_6*Omega^2*r2 + C_2_8*Omega^2*r2 + C_1_11*Omega^3*kon1*mu_9 + C_9_11*Omega^3*kon1*mu_1)/Omega^2;
xdot(33) = (Omega*koff1*mu_4 + 2*C_3_4*Omega^2*koff1 + C_3_9*Omega^2*kon1 - 2*C_3_3*Omega^3*kon1*mu_9 - 2*C_3_9*Omega^3*kon1*mu_3 + Omega^2*kon1*mu_3*mu_9)/Omega^2;
xdot(34) = -(Omega*koff1*mu_4 - kon1*(C_3_4*Omega^3*mu_9 - Omega*mu_4*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_4_9*Omega^2 + Omega^2*mu_4*mu_9) - Omega*mu_9*(C_3_4*Omega^2 + Omega^2*mu_3*mu_4) + C_3_9*Omega^3*mu_4 + C_4_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_4*mu_9) + C_3_4*Omega^2*koff1 - C_4_4*Omega^2*koff1 + C_3_9*Omega^2*kon1 - C_3_3*Omega^3*kon1*mu_9 - C_3_9*Omega^3*kon1*mu_3 + C_3_4*Omega^3*kon1*mu_9 + C_4_9*Omega^3*kon1*mu_3 + Omega^2*kon1*mu_3*mu_9)/Omega^2;
xdot(35) = (kon1*(C_3_5*Omega^3*mu_9 - Omega*mu_5*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) - Omega*mu_9*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_9*Omega^3*mu_5 + C_5_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_9) + kon2*(C_3_5*Omega^3*mu_10 - Omega*mu_5*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_10) + C_3_6*Omega^2*koff2 + C_4_5*Omega^2*koff1 - C_3_5*Omega^3*kon1*mu_9 - C_3_5*Omega^3*kon2*mu_10 - C_5_9*Omega^3*kon1*mu_3 - C_3_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(36) = (kon1*(C_3_6*Omega^3*mu_9 - Omega*mu_6*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_6_9*Omega^2 + Omega^2*mu_6*mu_9) - Omega*mu_9*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) + C_3_9*Omega^3*mu_6 + C_6_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_6*mu_9) - kon2*(C_3_5*Omega^3*mu_10 - Omega*mu_5*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_10) - C_3_6*Omega^2*koff2 + C_4_6*Omega^2*koff1 - C_3_6*Omega^3*kon1*mu_9 + C_3_5*Omega^3*kon2*mu_10 - C_6_9*Omega^3*kon1*mu_3 + C_3_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(37) = (kon1*(C_3_7*Omega^3*mu_9 - Omega*mu_7*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) - Omega*mu_9*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) + C_3_9*Omega^3*mu_7 + C_7_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_7*mu_9) + kon2*(C_3_7*Omega^3*mu_10 - Omega*mu_7*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) + C_3_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_7*mu_10) + C_3_8*Omega^2*koff2 + C_4_7*Omega^2*koff1 - C_3_7*Omega^3*kon1*mu_9 - C_3_7*Omega^3*kon2*mu_10 - C_7_9*Omega^3*kon1*mu_3 - C_3_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(38) = (kon1*(C_3_8*Omega^3*mu_9 - Omega*mu_8*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_8_9*Omega^2 + Omega^2*mu_8*mu_9) - Omega*mu_9*(C_3_8*Omega^2 + Omega^2*mu_3*mu_8) + C_3_9*Omega^3*mu_8 + C_8_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_8*mu_9) - kon2*(C_3_7*Omega^3*mu_10 - Omega*mu_7*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) + C_3_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_7*mu_10) - C_3_8*Omega^2*koff2 + C_4_8*Omega^2*koff1 - C_3_8*Omega^3*kon1*mu_9 + C_3_7*Omega^3*kon2*mu_10 - C_8_9*Omega^3*kon1*mu_3 + C_3_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(39) = (kon1*(C_1_3*Omega^3*mu_9 - Omega*mu_3*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_9*(C_1_3*Omega^2 + Omega^2*mu_1*mu_3) + C_1_9*Omega^3*mu_3 + C_3_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_3*mu_9) + Omega*koff1*mu_4 - C_3_9*Omega^2*betam1 + C_2_3*Omega^2*koff1 + C_3_4*Omega^2*koff1 + C_4_9*Omega^2*koff1 + C_3_9*Omega^2*kon1 - C_1_3*Omega^3*kon1*mu_9 - C_3_9*Omega^3*kon1*mu_1 - C_3_3*Omega^3*kon1*mu_9 - C_3_9*Omega^3*kon1*mu_3 - C_3_9*Omega^3*kon1*mu_9 - C_9_9*Omega^3*kon1*mu_3 + Omega^2*kon1*mu_3*mu_9)/Omega^2;
xdot(40) = -(C_3_10*Omega^2*betam2 - kon2*(C_3_7*Omega^3*mu_10 - Omega*mu_7*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) + C_3_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_7*mu_10) - kon1*(C_3_9*Omega^3*mu_10 - Omega*mu_9*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) + C_3_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_9*mu_10) - kon2*(C_3_5*Omega^3*mu_10 - Omega*mu_5*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_10) - C_3_6*Omega^2*koff2 - C_3_8*Omega^2*koff2 - C_4_10*Omega^2*koff1 + C_3_5*Omega^3*kon2*mu_10 + C_3_7*Omega^3*kon2*mu_10 + C_3_10*Omega^3*kon2*mu_5 + C_3_10*Omega^3*kon2*mu_7 + C_3_10*Omega^3*kon1*mu_9 + C_9_10*Omega^3*kon1*mu_3)/Omega^2;
xdot(41) = (kon1*(C_3_9*Omega^3*mu_11 - Omega*mu_9*(C_3_11*Omega^2 + Omega^2*mu_3*mu_11) - Omega*mu_3*(C_9_11*Omega^2 + Omega^2*mu_9*mu_11) - Omega*mu_11*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) + C_3_11*Omega^3*mu_9 + C_9_11*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_9*mu_11) + C_4_11*Omega^2*koff1 - C_3_11*Omega^2*myalpha + C_2_3*Omega^2*r1 + C_3_4*Omega^2*r1 + C_3_6*Omega^2*r2 + C_3_8*Omega^2*r2 - C_3_11*Omega^3*kon1*mu_9 - C_9_11*Omega^3*kon1*mu_3)/Omega^2;
xdot(42) = (Omega*koff1*mu_4 - 2*kon1*(C_3_4*Omega^3*mu_9 - Omega*mu_4*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_4_9*Omega^2 + Omega^2*mu_4*mu_9) - Omega*mu_9*(C_3_4*Omega^2 + Omega^2*mu_3*mu_4) + C_3_9*Omega^3*mu_4 + C_4_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_4*mu_9) - 2*C_4_4*Omega^2*koff1 + C_3_9*Omega^2*kon1 + 2*C_3_4*Omega^3*kon1*mu_9 + 2*C_4_9*Omega^3*kon1*mu_3 + Omega^2*kon1*mu_3*mu_9)/Omega^2;
xdot(43) = -(kon1*(C_3_5*Omega^3*mu_9 - Omega*mu_5*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) - Omega*mu_9*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_9*Omega^3*mu_5 + C_5_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_9) - kon2*(C_4_5*Omega^3*mu_10 - Omega*mu_5*(C_4_10*Omega^2 + Omega^2*mu_4*mu_10) - Omega*mu_4*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_4_5*Omega^2 + Omega^2*mu_4*mu_5) + C_4_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_4 + 3*Omega^3*mu_4*mu_5*mu_10) + C_4_5*Omega^2*koff1 - C_4_6*Omega^2*koff2 - C_3_5*Omega^3*kon1*mu_9 + C_4_5*Omega^3*kon2*mu_10 - C_5_9*Omega^3*kon1*mu_3 + C_4_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(44) = -(kon1*(C_3_6*Omega^3*mu_9 - Omega*mu_6*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_6_9*Omega^2 + Omega^2*mu_6*mu_9) - Omega*mu_9*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) + C_3_9*Omega^3*mu_6 + C_6_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_6*mu_9) + kon2*(C_4_5*Omega^3*mu_10 - Omega*mu_5*(C_4_10*Omega^2 + Omega^2*mu_4*mu_10) - Omega*mu_4*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_4_5*Omega^2 + Omega^2*mu_4*mu_5) + C_4_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_4 + 3*Omega^3*mu_4*mu_5*mu_10) + C_4_6*Omega^2*koff1 + C_4_6*Omega^2*koff2 - C_3_6*Omega^3*kon1*mu_9 - C_4_5*Omega^3*kon2*mu_10 - C_6_9*Omega^3*kon1*mu_3 - C_4_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(45) = -(kon1*(C_3_7*Omega^3*mu_9 - Omega*mu_7*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) - Omega*mu_9*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) + C_3_9*Omega^3*mu_7 + C_7_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_7*mu_9) - kon2*(C_4_7*Omega^3*mu_10 - Omega*mu_7*(C_4_10*Omega^2 + Omega^2*mu_4*mu_10) - Omega*mu_4*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_4_7*Omega^2 + Omega^2*mu_4*mu_7) + C_4_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_4 + 3*Omega^3*mu_4*mu_7*mu_10) + C_4_7*Omega^2*koff1 - C_4_8*Omega^2*koff2 - C_3_7*Omega^3*kon1*mu_9 + C_4_7*Omega^3*kon2*mu_10 - C_7_9*Omega^3*kon1*mu_3 + C_4_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(46) = -(kon1*(C_3_8*Omega^3*mu_9 - Omega*mu_8*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_8_9*Omega^2 + Omega^2*mu_8*mu_9) - Omega*mu_9*(C_3_8*Omega^2 + Omega^2*mu_3*mu_8) + C_3_9*Omega^3*mu_8 + C_8_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_8*mu_9) + kon2*(C_4_7*Omega^3*mu_10 - Omega*mu_7*(C_4_10*Omega^2 + Omega^2*mu_4*mu_10) - Omega*mu_4*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_4_7*Omega^2 + Omega^2*mu_4*mu_7) + C_4_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_4 + 3*Omega^3*mu_4*mu_7*mu_10) + C_4_8*Omega^2*koff1 + C_4_8*Omega^2*koff2 - C_3_8*Omega^3*kon1*mu_9 - C_4_7*Omega^3*kon2*mu_10 - C_8_9*Omega^3*kon1*mu_3 - C_4_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(47) = -(Omega*koff1*mu_4 - kon1*(C_3_4*Omega^3*mu_9 - Omega*mu_4*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_4_9*Omega^2 + Omega^2*mu_4*mu_9) - Omega*mu_9*(C_3_4*Omega^2 + Omega^2*mu_3*mu_4) + C_3_9*Omega^3*mu_4 + C_4_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_4*mu_9) - kon1*(C_1_4*Omega^3*mu_9 - Omega*mu_4*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_4_9*Omega^2 + Omega^2*mu_4*mu_9) - Omega*mu_9*(C_1_4*Omega^2 + Omega^2*mu_1*mu_4) + C_1_9*Omega^3*mu_4 + C_4_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_4*mu_9) + C_4_9*Omega^2*betam1 - C_2_4*Omega^2*koff1 - C_4_4*Omega^2*koff1 + C_4_9*Omega^2*koff1 + C_3_9*Omega^2*kon1 + C_1_4*Omega^3*kon1*mu_9 + C_3_4*Omega^3*kon1*mu_9 - C_3_9*Omega^3*kon1*mu_9 + C_4_9*Omega^3*kon1*mu_1 + C_4_9*Omega^3*kon1*mu_3 - C_9_9*Omega^3*kon1*mu_3 + Omega^2*kon1*mu_3*mu_9)/Omega^2;
xdot(48) = -(kon1*(C_3_9*Omega^3*mu_10 - Omega*mu_9*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) + C_3_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_9*mu_10) - kon2*(C_4_7*Omega^3*mu_10 - Omega*mu_7*(C_4_10*Omega^2 + Omega^2*mu_4*mu_10) - Omega*mu_4*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_4_7*Omega^2 + Omega^2*mu_4*mu_7) + C_4_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_4 + 3*Omega^3*mu_4*mu_7*mu_10) - kon2*(C_4_5*Omega^3*mu_10 - Omega*mu_5*(C_4_10*Omega^2 + Omega^2*mu_4*mu_10) - Omega*mu_4*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_10*(C_4_5*Omega^2 + Omega^2*mu_4*mu_5) + C_4_10*Omega^3*mu_5 + C_5_10*Omega^3*mu_4 + 3*Omega^3*mu_4*mu_5*mu_10) + C_4_10*Omega^2*betam2 - C_4_6*Omega^2*koff2 - C_4_8*Omega^2*koff2 + C_4_10*Omega^2*koff1 + C_4_5*Omega^3*kon2*mu_10 + C_4_7*Omega^3*kon2*mu_10 - C_3_10*Omega^3*kon1*mu_9 + C_4_10*Omega^3*kon2*mu_5 + C_4_10*Omega^3*kon2*mu_7 - C_9_10*Omega^3*kon1*mu_3)/Omega^2;
xdot(49) = (C_2_4*Omega^2*r1 - C_4_11*Omega^2*koff1 - C_4_11*Omega^2*myalpha - kon1*(C_3_9*Omega^3*mu_11 - Omega*mu_9*(C_3_11*Omega^2 + Omega^2*mu_3*mu_11) - Omega*mu_3*(C_9_11*Omega^2 + Omega^2*mu_9*mu_11) - Omega*mu_11*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) + C_3_11*Omega^3*mu_9 + C_9_11*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_9*mu_11) + C_4_4*Omega^2*r1 + C_4_6*Omega^2*r2 + C_4_8*Omega^2*r2 + C_3_11*Omega^3*kon1*mu_9 + C_9_11*Omega^3*kon1*mu_3)/Omega^2;
xdot(50) = (Omega*koff2*mu_6 + 2*C_5_6*Omega^2*koff2 + C_5_10*Omega^2*kon2 - 2*C_5_5*Omega^3*kon2*mu_10 - 2*C_5_10*Omega^3*kon2*mu_5 + Omega^2*kon2*mu_5*mu_10)/Omega^2;
xdot(51) = -(Omega*koff2*mu_6 - kon2*(C_5_6*Omega^3*mu_10 - Omega*mu_6*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_6_10*Omega^2 + Omega^2*mu_6*mu_10) - Omega*mu_10*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) + C_5_10*Omega^3*mu_6 + C_6_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_6*mu_10) + C_5_6*Omega^2*koff2 - C_6_6*Omega^2*koff2 + C_5_10*Omega^2*kon2 - C_5_5*Omega^3*kon2*mu_10 + C_5_6*Omega^3*kon2*mu_10 - C_5_10*Omega^3*kon2*mu_5 + C_6_10*Omega^3*kon2*mu_5 + Omega^2*kon2*mu_5*mu_10)/Omega^2;
xdot(52) = (2*kon2*(C_5_7*Omega^3*mu_10 - Omega*mu_7*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) + C_5_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_7*mu_10) + C_5_8*Omega^2*koff2 + C_6_7*Omega^2*koff2 - 2*C_5_7*Omega^3*kon2*mu_10 - C_5_10*Omega^3*kon2*mu_7 - C_7_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(53) = -(kon2*(C_5_7*Omega^3*mu_10 - Omega*mu_7*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) + C_5_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_7*mu_10) - kon2*(C_5_8*Omega^3*mu_10 - Omega*mu_8*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_8_10*Omega^2 + Omega^2*mu_8*mu_10) - Omega*mu_10*(C_5_8*Omega^2 + Omega^2*mu_5*mu_8) + C_5_10*Omega^3*mu_8 + C_8_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_8*mu_10) + C_5_8*Omega^2*koff2 - C_6_8*Omega^2*koff2 - C_5_7*Omega^3*kon2*mu_10 + C_5_8*Omega^3*kon2*mu_10 - C_5_10*Omega^3*kon2*mu_7 + C_8_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(54) = -(C_5_9*Omega^2*betam1 - kon1*(C_3_5*Omega^3*mu_9 - Omega*mu_5*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) - Omega*mu_9*(C_3_5*Omega^2 + Omega^2*mu_3*mu_5) + C_3_9*Omega^3*mu_5 + C_5_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_5*mu_9) - kon2*(C_5_9*Omega^3*mu_10 - Omega*mu_9*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) + C_5_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_9*mu_10) - kon1*(C_1_5*Omega^3*mu_9 - Omega*mu_5*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) - Omega*mu_9*(C_1_5*Omega^2 + Omega^2*mu_1*mu_5) + C_1_9*Omega^3*mu_5 + C_5_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_5*mu_9) - C_2_5*Omega^2*koff1 - C_4_5*Omega^2*koff1 - C_6_9*Omega^2*koff2 + C_1_5*Omega^3*kon1*mu_9 + C_3_5*Omega^3*kon1*mu_9 + C_5_9*Omega^3*kon1*mu_1 + C_5_9*Omega^3*kon1*mu_3 + C_5_9*Omega^3*kon2*mu_10 + C_9_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(55) = (kon2*(C_5_7*Omega^3*mu_10 - Omega*mu_7*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) + C_5_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_7*mu_10) + Omega*koff2*mu_6 - C_5_10*Omega^2*betam2 + C_5_6*Omega^2*koff2 + C_5_8*Omega^2*koff2 + C_6_10*Omega^2*koff2 + C_5_10*Omega^2*kon2 - C_5_5*Omega^3*kon2*mu_10 - C_5_7*Omega^3*kon2*mu_10 - C_5_10*Omega^3*kon2*mu_5 - C_5_10*Omega^3*kon2*mu_7 - C_5_10*Omega^3*kon2*mu_10 - C_10_10*Omega^3*kon2*mu_5 + Omega^2*kon2*mu_5*mu_10)/Omega^2;
xdot(56) = (kon2*(C_5_10*Omega^3*mu_11 - Omega*mu_10*(C_5_11*Omega^2 + Omega^2*mu_5*mu_11) - Omega*mu_5*(C_10_11*Omega^2 + Omega^2*mu_10*mu_11) - Omega*mu_11*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) + C_5_11*Omega^3*mu_10 + C_10_11*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_10*mu_11) + C_6_11*Omega^2*koff2 - C_5_11*Omega^2*myalpha + C_2_5*Omega^2*r1 + C_4_5*Omega^2*r1 + C_5_6*Omega^2*r2 + C_5_8*Omega^2*r2 - C_5_11*Omega^3*kon2*mu_10 - C_10_11*Omega^3*kon2*mu_5)/Omega^2;
xdot(57) = (Omega*koff2*mu_6 - 2*kon2*(C_5_6*Omega^3*mu_10 - Omega*mu_6*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_6_10*Omega^2 + Omega^2*mu_6*mu_10) - Omega*mu_10*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) + C_5_10*Omega^3*mu_6 + C_6_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_6*mu_10) - 2*C_6_6*Omega^2*koff2 + C_5_10*Omega^2*kon2 + 2*C_5_6*Omega^3*kon2*mu_10 + 2*C_6_10*Omega^3*kon2*mu_5 + Omega^2*kon2*mu_5*mu_10)/Omega^2;
xdot(58) = -(kon2*(C_5_7*Omega^3*mu_10 - Omega*mu_7*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) + C_5_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_7*mu_10) - kon2*(C_6_7*Omega^3*mu_10 - Omega*mu_7*(C_6_10*Omega^2 + Omega^2*mu_6*mu_10) - Omega*mu_6*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_6_7*Omega^2 + Omega^2*mu_6*mu_7) + C_6_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_6 + 3*Omega^3*mu_6*mu_7*mu_10) + C_6_7*Omega^2*koff2 - C_6_8*Omega^2*koff2 - C_5_7*Omega^3*kon2*mu_10 + C_6_7*Omega^3*kon2*mu_10 + C_6_10*Omega^3*kon2*mu_7 - C_7_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(59) = (C_5_8*Omega^3*kon2*mu_10 - kon2*(C_6_7*Omega^3*mu_10 - Omega*mu_7*(C_6_10*Omega^2 + Omega^2*mu_6*mu_10) - Omega*mu_6*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_6_7*Omega^2 + Omega^2*mu_6*mu_7) + C_6_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_6 + 3*Omega^3*mu_6*mu_7*mu_10) - 2*C_6_8*Omega^2*koff2 - kon2*(C_5_8*Omega^3*mu_10 - Omega*mu_8*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_8_10*Omega^2 + Omega^2*mu_8*mu_10) - Omega*mu_10*(C_5_8*Omega^2 + Omega^2*mu_5*mu_8) + C_5_10*Omega^3*mu_8 + C_8_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_8*mu_10) + C_6_7*Omega^3*kon2*mu_10 + C_6_10*Omega^3*kon2*mu_7 + C_8_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(60) = -(kon2*(C_5_9*Omega^3*mu_10 - Omega*mu_9*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) + C_5_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_9*mu_10) - kon1*(C_3_6*Omega^3*mu_9 - Omega*mu_6*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_6_9*Omega^2 + Omega^2*mu_6*mu_9) - Omega*mu_9*(C_3_6*Omega^2 + Omega^2*mu_3*mu_6) + C_3_9*Omega^3*mu_6 + C_6_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_6*mu_9) - kon1*(C_1_6*Omega^3*mu_9 - Omega*mu_6*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_6_9*Omega^2 + Omega^2*mu_6*mu_9) - Omega*mu_9*(C_1_6*Omega^2 + Omega^2*mu_1*mu_6) + C_1_9*Omega^3*mu_6 + C_6_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_6*mu_9) + C_6_9*Omega^2*betam1 - C_2_6*Omega^2*koff1 - C_4_6*Omega^2*koff1 + C_6_9*Omega^2*koff2 + C_1_6*Omega^3*kon1*mu_9 + C_3_6*Omega^3*kon1*mu_9 - C_5_9*Omega^3*kon2*mu_10 + C_6_9*Omega^3*kon1*mu_1 + C_6_9*Omega^3*kon1*mu_3 - C_9_10*Omega^3*kon2*mu_5)/Omega^2;
xdot(61) = -(Omega*koff2*mu_6 - kon2*(C_6_7*Omega^3*mu_10 - Omega*mu_7*(C_6_10*Omega^2 + Omega^2*mu_6*mu_10) - Omega*mu_6*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_6_7*Omega^2 + Omega^2*mu_6*mu_7) + C_6_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_6 + 3*Omega^3*mu_6*mu_7*mu_10) - kon2*(C_5_6*Omega^3*mu_10 - Omega*mu_6*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_6_10*Omega^2 + Omega^2*mu_6*mu_10) - Omega*mu_10*(C_5_6*Omega^2 + Omega^2*mu_5*mu_6) + C_5_10*Omega^3*mu_6 + C_6_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_6*mu_10) + C_6_10*Omega^2*betam2 - C_6_6*Omega^2*koff2 - C_6_8*Omega^2*koff2 + C_6_10*Omega^2*koff2 + C_5_10*Omega^2*kon2 + C_5_6*Omega^3*kon2*mu_10 + C_6_7*Omega^3*kon2*mu_10 - C_5_10*Omega^3*kon2*mu_10 + C_6_10*Omega^3*kon2*mu_5 + C_6_10*Omega^3*kon2*mu_7 - C_10_10*Omega^3*kon2*mu_5 + Omega^2*kon2*mu_5*mu_10)/Omega^2;
xdot(62) = (C_2_6*Omega^2*r1 - C_6_11*Omega^2*koff2 - C_6_11*Omega^2*myalpha - kon2*(C_5_10*Omega^3*mu_11 - Omega*mu_10*(C_5_11*Omega^2 + Omega^2*mu_5*mu_11) - Omega*mu_5*(C_10_11*Omega^2 + Omega^2*mu_10*mu_11) - Omega*mu_11*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) + C_5_11*Omega^3*mu_10 + C_10_11*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_10*mu_11) + C_4_6*Omega^2*r1 + C_6_6*Omega^2*r2 + C_6_8*Omega^2*r2 + C_5_11*Omega^3*kon2*mu_10 + C_10_11*Omega^3*kon2*mu_5)/Omega^2;
xdot(63) = (Omega*koff2*mu_8 + 2*C_7_8*Omega^2*koff2 + C_7_10*Omega^2*kon2 - 2*C_7_7*Omega^3*kon2*mu_10 - 2*C_7_10*Omega^3*kon2*mu_7 + Omega^2*kon2*mu_7*mu_10)/Omega^2;
xdot(64) = -(Omega*koff2*mu_8 - kon2*(C_7_8*Omega^3*mu_10 - Omega*mu_8*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_7*(C_8_10*Omega^2 + Omega^2*mu_8*mu_10) - Omega*mu_10*(C_7_8*Omega^2 + Omega^2*mu_7*mu_8) + C_7_10*Omega^3*mu_8 + C_8_10*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_8*mu_10) + C_7_8*Omega^2*koff2 - C_8_8*Omega^2*koff2 + C_7_10*Omega^2*kon2 - C_7_7*Omega^3*kon2*mu_10 + C_7_8*Omega^3*kon2*mu_10 - C_7_10*Omega^3*kon2*mu_7 + C_8_10*Omega^3*kon2*mu_7 + Omega^2*kon2*mu_7*mu_10)/Omega^2;
xdot(65) = -(C_7_9*Omega^2*betam1 - kon1*(C_3_7*Omega^3*mu_9 - Omega*mu_7*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) - Omega*mu_9*(C_3_7*Omega^2 + Omega^2*mu_3*mu_7) + C_3_9*Omega^3*mu_7 + C_7_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_7*mu_9) - kon2*(C_7_9*Omega^3*mu_10 - Omega*mu_9*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_7*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) + C_7_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_9*mu_10) - kon1*(C_1_7*Omega^3*mu_9 - Omega*mu_7*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) - Omega*mu_9*(C_1_7*Omega^2 + Omega^2*mu_1*mu_7) + C_1_9*Omega^3*mu_7 + C_7_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_7*mu_9) - C_2_7*Omega^2*koff1 - C_4_7*Omega^2*koff1 - C_8_9*Omega^2*koff2 + C_1_7*Omega^3*kon1*mu_9 + C_3_7*Omega^3*kon1*mu_9 + C_7_9*Omega^3*kon1*mu_1 + C_7_9*Omega^3*kon1*mu_3 + C_7_9*Omega^3*kon2*mu_10 + C_9_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(66) = (kon2*(C_5_7*Omega^3*mu_10 - Omega*mu_7*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_10*(C_5_7*Omega^2 + Omega^2*mu_5*mu_7) + C_5_10*Omega^3*mu_7 + C_7_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_7*mu_10) + Omega*koff2*mu_8 - C_7_10*Omega^2*betam2 + C_6_7*Omega^2*koff2 + C_7_8*Omega^2*koff2 + C_8_10*Omega^2*koff2 + C_7_10*Omega^2*kon2 - C_5_7*Omega^3*kon2*mu_10 - C_7_7*Omega^3*kon2*mu_10 - C_7_10*Omega^3*kon2*mu_5 - C_7_10*Omega^3*kon2*mu_7 - C_7_10*Omega^3*kon2*mu_10 - C_10_10*Omega^3*kon2*mu_7 + Omega^2*kon2*mu_7*mu_10)/Omega^2;
xdot(67) = (kon2*(C_7_10*Omega^3*mu_11 - Omega*mu_10*(C_7_11*Omega^2 + Omega^2*mu_7*mu_11) - Omega*mu_7*(C_10_11*Omega^2 + Omega^2*mu_10*mu_11) - Omega*mu_11*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) + C_7_11*Omega^3*mu_10 + C_10_11*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_10*mu_11) + C_8_11*Omega^2*koff2 - C_7_11*Omega^2*myalpha + C_2_7*Omega^2*r1 + C_4_7*Omega^2*r1 + C_6_7*Omega^2*r2 + C_7_8*Omega^2*r2 - C_7_11*Omega^3*kon2*mu_10 - C_10_11*Omega^3*kon2*mu_7)/Omega^2;
xdot(68) = (Omega*koff2*mu_8 - 2*kon2*(C_7_8*Omega^3*mu_10 - Omega*mu_8*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_7*(C_8_10*Omega^2 + Omega^2*mu_8*mu_10) - Omega*mu_10*(C_7_8*Omega^2 + Omega^2*mu_7*mu_8) + C_7_10*Omega^3*mu_8 + C_8_10*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_8*mu_10) - 2*C_8_8*Omega^2*koff2 + C_7_10*Omega^2*kon2 + 2*C_7_8*Omega^3*kon2*mu_10 + 2*C_8_10*Omega^3*kon2*mu_7 + Omega^2*kon2*mu_7*mu_10)/Omega^2;
xdot(69) = -(kon2*(C_7_9*Omega^3*mu_10 - Omega*mu_9*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_7*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) + C_7_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_9*mu_10) - kon1*(C_3_8*Omega^3*mu_9 - Omega*mu_8*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) - Omega*mu_3*(C_8_9*Omega^2 + Omega^2*mu_8*mu_9) - Omega*mu_9*(C_3_8*Omega^2 + Omega^2*mu_3*mu_8) + C_3_9*Omega^3*mu_8 + C_8_9*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_8*mu_9) - kon1*(C_1_8*Omega^3*mu_9 - Omega*mu_8*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) - Omega*mu_1*(C_8_9*Omega^2 + Omega^2*mu_8*mu_9) - Omega*mu_9*(C_1_8*Omega^2 + Omega^2*mu_1*mu_8) + C_1_9*Omega^3*mu_8 + C_8_9*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_8*mu_9) + C_8_9*Omega^2*betam1 - C_2_8*Omega^2*koff1 - C_4_8*Omega^2*koff1 + C_8_9*Omega^2*koff2 + C_1_8*Omega^3*kon1*mu_9 + C_3_8*Omega^3*kon1*mu_9 - C_7_9*Omega^3*kon2*mu_10 + C_8_9*Omega^3*kon1*mu_1 + C_8_9*Omega^3*kon1*mu_3 - C_9_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(70) = -(Omega*koff2*mu_8 - kon2*(C_7_8*Omega^3*mu_10 - Omega*mu_8*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_7*(C_8_10*Omega^2 + Omega^2*mu_8*mu_10) - Omega*mu_10*(C_7_8*Omega^2 + Omega^2*mu_7*mu_8) + C_7_10*Omega^3*mu_8 + C_8_10*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_8*mu_10) - kon2*(C_5_8*Omega^3*mu_10 - Omega*mu_8*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_8_10*Omega^2 + Omega^2*mu_8*mu_10) - Omega*mu_10*(C_5_8*Omega^2 + Omega^2*mu_5*mu_8) + C_5_10*Omega^3*mu_8 + C_8_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_8*mu_10) + C_8_10*Omega^2*betam2 - C_6_8*Omega^2*koff2 - C_8_8*Omega^2*koff2 + C_8_10*Omega^2*koff2 + C_7_10*Omega^2*kon2 + C_5_8*Omega^3*kon2*mu_10 + C_7_8*Omega^3*kon2*mu_10 - C_7_10*Omega^3*kon2*mu_10 + C_8_10*Omega^3*kon2*mu_5 + C_8_10*Omega^3*kon2*mu_7 - C_10_10*Omega^3*kon2*mu_7 + Omega^2*kon2*mu_7*mu_10)/Omega^2;
xdot(71) = (C_2_8*Omega^2*r1 - C_8_11*Omega^2*koff2 - C_8_11*Omega^2*myalpha - kon2*(C_7_10*Omega^3*mu_11 - Omega*mu_10*(C_7_11*Omega^2 + Omega^2*mu_7*mu_11) - Omega*mu_7*(C_10_11*Omega^2 + Omega^2*mu_10*mu_11) - Omega*mu_11*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) + C_7_11*Omega^3*mu_10 + C_10_11*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_10*mu_11) + C_4_8*Omega^2*r1 + C_6_8*Omega^2*r2 + C_8_8*Omega^2*r2 + C_7_11*Omega^3*kon2*mu_10 + C_10_11*Omega^3*kon2*mu_7)/Omega^2;
xdot(72) = (16*beta1 + Omega*betam1*mu_9 + Omega*koff1*mu_2 + Omega*koff1*mu_4 - 2*C_9_9*Omega^2*betam1 + 2*C_2_9*Omega^2*koff1 + 2*C_4_9*Omega^2*koff1 + C_1_9*Omega^2*kon1 + C_3_9*Omega^2*kon1 - 2*C_1_9*Omega^3*kon1*mu_9 - 2*C_3_9*Omega^3*kon1*mu_9 - 2*C_9_9*Omega^3*kon1*mu_1 - 2*C_9_9*Omega^3*kon1*mu_3 + Omega^2*kon1*mu_1*mu_9 + Omega^2*kon1*mu_3*mu_9)/Omega^2;
xdot(73) = -(C_9_10*Omega^2*betam1 - kon1*(C_3_9*Omega^3*mu_10 - Omega*mu_9*(C_3_10*Omega^2 + Omega^2*mu_3*mu_10) - Omega*mu_3*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) + C_3_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_9*mu_10) - kon2*(C_5_9*Omega^3*mu_10 - Omega*mu_9*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) - Omega*mu_5*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_5_9*Omega^2 + Omega^2*mu_5*mu_9) + C_5_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_9*mu_10) - kon2*(C_7_9*Omega^3*mu_10 - Omega*mu_9*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) - Omega*mu_7*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_7_9*Omega^2 + Omega^2*mu_7*mu_9) + C_7_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_9*mu_10) - kon1*(C_1_9*Omega^3*mu_10 - Omega*mu_9*(C_1_10*Omega^2 + Omega^2*mu_1*mu_10) - Omega*mu_1*(C_9_10*Omega^2 + Omega^2*mu_9*mu_10) - Omega*mu_10*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) + C_1_10*Omega^3*mu_9 + C_9_10*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_9*mu_10) + C_9_10*Omega^2*betam2 - C_6_9*Omega^2*koff2 - C_8_9*Omega^2*koff2 - C_2_10*Omega^2*koff1 - C_4_10*Omega^2*koff1 + C_5_9*Omega^3*kon2*mu_10 + C_7_9*Omega^3*kon2*mu_10 + C_1_10*Omega^3*kon1*mu_9 + C_3_10*Omega^3*kon1*mu_9 + C_9_10*Omega^3*kon1*mu_1 + C_9_10*Omega^3*kon1*mu_3 + C_9_10*Omega^3*kon2*mu_5 + C_9_10*Omega^3*kon2*mu_7)/Omega^2;
xdot(74) = (kon1*(C_1_9*Omega^3*mu_11 - Omega*mu_9*(C_1_11*Omega^2 + Omega^2*mu_1*mu_11) - Omega*mu_1*(C_9_11*Omega^2 + Omega^2*mu_9*mu_11) - Omega*mu_11*(C_1_9*Omega^2 + Omega^2*mu_1*mu_9) + C_1_11*Omega^3*mu_9 + C_9_11*Omega^3*mu_1 + 3*Omega^3*mu_1*mu_9*mu_11) + kon1*(C_3_9*Omega^3*mu_11 - Omega*mu_9*(C_3_11*Omega^2 + Omega^2*mu_3*mu_11) - Omega*mu_3*(C_9_11*Omega^2 + Omega^2*mu_9*mu_11) - Omega*mu_11*(C_3_9*Omega^2 + Omega^2*mu_3*mu_9) + C_3_11*Omega^3*mu_9 + C_9_11*Omega^3*mu_3 + 3*Omega^3*mu_3*mu_9*mu_11) - C_9_11*Omega^2*betam1 + C_2_11*Omega^2*koff1 + C_4_11*Omega^2*koff1 - C_9_11*Omega^2*myalpha + C_2_9*Omega^2*r1 + C_4_9*Omega^2*r1 + C_6_9*Omega^2*r2 + C_8_9*Omega^2*r2 - C_1_11*Omega^3*kon1*mu_9 - C_3_11*Omega^3*kon1*mu_9 - C_9_11*Omega^3*kon1*mu_1 - C_9_11*Omega^3*kon1*mu_3)/Omega^2;
xdot(75) = (144*beta2 + Omega*betam2*mu_10 + Omega*koff2*mu_6 + Omega*koff2*mu_8 - 2*C_10_10*Omega^2*betam2 + 2*C_6_10*Omega^2*koff2 + 2*C_8_10*Omega^2*koff2 + C_5_10*Omega^2*kon2 + C_7_10*Omega^2*kon2 - 2*C_5_10*Omega^3*kon2*mu_10 - 2*C_7_10*Omega^3*kon2*mu_10 - 2*C_10_10*Omega^3*kon2*mu_5 - 2*C_10_10*Omega^3*kon2*mu_7 + Omega^2*kon2*mu_5*mu_10 + Omega^2*kon2*mu_7*mu_10)/Omega^2;
xdot(76) = (kon2*(C_5_10*Omega^3*mu_11 - Omega*mu_10*(C_5_11*Omega^2 + Omega^2*mu_5*mu_11) - Omega*mu_5*(C_10_11*Omega^2 + Omega^2*mu_10*mu_11) - Omega*mu_11*(C_5_10*Omega^2 + Omega^2*mu_5*mu_10) + C_5_11*Omega^3*mu_10 + C_10_11*Omega^3*mu_5 + 3*Omega^3*mu_5*mu_10*mu_11) + kon2*(C_7_10*Omega^3*mu_11 - Omega*mu_10*(C_7_11*Omega^2 + Omega^2*mu_7*mu_11) - Omega*mu_7*(C_10_11*Omega^2 + Omega^2*mu_10*mu_11) - Omega*mu_11*(C_7_10*Omega^2 + Omega^2*mu_7*mu_10) + C_7_11*Omega^3*mu_10 + C_10_11*Omega^3*mu_7 + 3*Omega^3*mu_7*mu_10*mu_11) - C_10_11*Omega^2*betam2 + C_6_11*Omega^2*koff2 + C_8_11*Omega^2*koff2 - C_10_11*Omega^2*myalpha + C_2_10*Omega^2*r1 + C_4_10*Omega^2*r1 + C_6_10*Omega^2*r2 + C_8_10*Omega^2*r2 - C_5_11*Omega^3*kon2*mu_10 - C_7_11*Omega^3*kon2*mu_10 - C_10_11*Omega^3*kon2*mu_5 - C_10_11*Omega^3*kon2*mu_7)/Omega^2;
xdot(77) = (Omega*mu_11*myalpha + Omega*mu_2*r1 + Omega*mu_4*r1 + Omega*mu_6*r2 + Omega*mu_8*r2 - 2*C_11_11*Omega^2*myalpha + 2*C_2_11*Omega^2*r1 + 2*C_4_11*Omega^2*r1 + 2*C_6_11*Omega^2*r2 + 2*C_8_11*Omega^2*r2)/Omega^2;
% INITIAL CONDITIONS

x0 = sym(zeros(size(x)));

x0(1) = (indmu1*kmu01 - fmu01*(indmu1 - 1))/Omega;
x0(2) = (indmu2*kmu02 - fmu02*(indmu2 - 1))/Omega;
x0(3) = (indmu3*kmu03 - fmu03*(indmu3 - 1))/Omega;
x0(4) = (indmu4*kmu04 - fmu04*(indmu4 - 1))/Omega;
x0(5) = (indmu5*kmu05 - fmu05*(indmu5 - 1))/Omega;
x0(6) = (indmu6*kmu06 - fmu06*(indmu6 - 1))/Omega;
x0(7) = (indmu7*kmu07 - fmu07*(indmu7 - 1))/Omega;
x0(8) = (indmu8*kmu08 - fmu08*(indmu8 - 1))/Omega;
x0(9) = (indmu9*kmu09 - fmu09*(indmu9 - 1))/Omega;
x0(10) = (indmu10*kmu010 - fmu010*(indmu10 - 1))/Omega;
x0(11) = (indmu11*kmu011 - fmu011*(indmu11 - 1))/Omega;
x0(12) = (indC1*kC01 - fC01*(indC1 - 1))/Omega^2;
x0(13) = (indC2*kC02 - fC02*(indC2 - 1))/Omega^2;
x0(14) = (indC3*kC03 - fC03*(indC3 - 1))/Omega^2;
x0(15) = (indC4*kC04 - fC04*(indC4 - 1))/Omega^2;
x0(16) = (indC5*kC05 - fC05*(indC5 - 1))/Omega^2;
x0(17) = (indC6*kC06 - fC06*(indC6 - 1))/Omega^2;
x0(18) = (indC7*kC07 - fC07*(indC7 - 1))/Omega^2;
x0(19) = (indC8*kC08 - fC08*(indC8 - 1))/Omega^2;
x0(20) = (indC9*kC09 - fC09*(indC9 - 1))/Omega^2;
x0(21) = (indC10*kC010 - fC010*(indC10 - 1))/Omega^2;
x0(22) = (indC11*kC011 - fC011*(indC11 - 1))/Omega^2;
x0(23) = (indC12*kC012 - fC012*(indC12 - 1))/Omega^2;
x0(24) = (indC13*kC013 - fC013*(indC13 - 1))/Omega^2;
x0(25) = (indC14*kC014 - fC014*(indC14 - 1))/Omega^2;
x0(26) = (indC15*kC015 - fC015*(indC15 - 1))/Omega^2;
x0(27) = (indC16*kC016 - fC016*(indC16 - 1))/Omega^2;
x0(28) = (indC17*kC017 - fC017*(indC17 - 1))/Omega^2;
x0(29) = (indC18*kC018 - fC018*(indC18 - 1))/Omega^2;
x0(30) = (indC19*kC019 - fC019*(indC19 - 1))/Omega^2;
x0(31) = (indC20*kC020 - fC020*(indC20 - 1))/Omega^2;
x0(32) = (indC21*kC021 - fC021*(indC21 - 1))/Omega^2;
x0(33) = (indC22*kC022 - fC022*(indC22 - 1))/Omega^2;
x0(34) = (indC23*kC023 - fC023*(indC23 - 1))/Omega^2;
x0(35) = (indC24*kC024 - fC024*(indC24 - 1))/Omega^2;
x0(36) = (indC25*kC025 - fC025*(indC25 - 1))/Omega^2;
x0(37) = (indC26*kC026 - fC026*(indC26 - 1))/Omega^2;
x0(38) = (indC27*kC027 - fC027*(indC27 - 1))/Omega^2;
x0(39) = (indC28*kC028 - fC028*(indC28 - 1))/Omega^2;
x0(40) = (indC29*kC029 - fC029*(indC29 - 1))/Omega^2;
x0(41) = (indC30*kC030 - fC030*(indC30 - 1))/Omega^2;
x0(42) = (indC31*kC031 - fC031*(indC31 - 1))/Omega^2;
x0(43) = (indC32*kC032 - fC032*(indC32 - 1))/Omega^2;
x0(44) = (indC33*kC033 - fC033*(indC33 - 1))/Omega^2;
x0(45) = (indC34*kC034 - fC034*(indC34 - 1))/Omega^2;
x0(46) = (indC35*kC035 - fC035*(indC35 - 1))/Omega^2;
x0(47) = (indC36*kC036 - fC036*(indC36 - 1))/Omega^2;
x0(48) = (indC37*kC037 - fC037*(indC37 - 1))/Omega^2;
x0(49) = (indC38*kC038 - fC038*(indC38 - 1))/Omega^2;
x0(50) = (indC39*kC039 - fC039*(indC39 - 1))/Omega^2;
x0(51) = (indC40*kC040 - fC040*(indC40 - 1))/Omega^2;
x0(52) = (indC41*kC041 - fC041*(indC41 - 1))/Omega^2;
x0(53) = (indC42*kC042 - fC042*(indC42 - 1))/Omega^2;
x0(54) = (indC43*kC043 - fC043*(indC43 - 1))/Omega^2;
x0(55) = (indC44*kC044 - fC044*(indC44 - 1))/Omega^2;
x0(56) = (indC45*kC045 - fC045*(indC45 - 1))/Omega^2;
x0(57) = (indC46*kC046 - fC046*(indC46 - 1))/Omega^2;
x0(58) = (indC47*kC047 - fC047*(indC47 - 1))/Omega^2;
x0(59) = (indC48*kC048 - fC048*(indC48 - 1))/Omega^2;
x0(60) = (indC49*kC049 - fC049*(indC49 - 1))/Omega^2;
x0(61) = (indC50*kC050 - fC050*(indC50 - 1))/Omega^2;
x0(62) = (indC51*kC051 - fC051*(indC51 - 1))/Omega^2;
x0(63) = (indC52*kC052 - fC052*(indC52 - 1))/Omega^2;
x0(64) = (indC53*kC053 - fC053*(indC53 - 1))/Omega^2;
x0(65) = (indC54*kC054 - fC054*(indC54 - 1))/Omega^2;
x0(66) = (indC55*kC055 - fC055*(indC55 - 1))/Omega^2;
x0(67) = (indC56*kC056 - fC056*(indC56 - 1))/Omega^2;
x0(68) = (indC57*kC057 - fC057*(indC57 - 1))/Omega^2;
x0(69) = (indC58*kC058 - fC058*(indC58 - 1))/Omega^2;
x0(70) = (indC59*kC059 - fC059*(indC59 - 1))/Omega^2;
x0(71) = (indC60*kC060 - fC060*(indC60 - 1))/Omega^2;
x0(72) = (indC61*kC061 - fC061*(indC61 - 1))/Omega^2;
x0(73) = (indC62*kC062 - fC062*(indC62 - 1))/Omega^2;
x0(74) = (indC63*kC063 - fC063*(indC63 - 1))/Omega^2;
x0(75) = (indC64*kC064 - fC064*(indC64 - 1))/Omega^2;
x0(76) = (indC65*kC065 - fC065*(indC65 - 1))/Omega^2;
x0(77) = (indC66*kC066 - fC066*(indC66 - 1))/Omega^2;

% OBSERVABLES

y = sym(zeros(2,1));

y(1) = mu_11;
y(2) = C_11_11;

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