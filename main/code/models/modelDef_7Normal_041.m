syms A0000 A0002 A0020 A0022 A0200 A0202 A0220 A0222 A2000 A2002 A2020 A2022 A2200 A2202 A2220 A2222 T2 kon2 koff2 r2 myalpha beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0000; A0002; A0020; A0022; A0200; A0202; A0220; A0222; A2000; A2002; A2020; A2022; A2200; A2202; A2220; A2222; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8) + x(9) + x(10) + x(11) + x(12) + x(13) + x(14) + x(15) + x(16)) == 1);
System.parameter.variable = [kon2; koff2; r2; myalpha; beta2; betam2];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A0002;
System.reaction(1).product = [A0000,T2];
System.reaction(1).propensity = koff2 * A0002;

System.reaction(2).educt = A0020;
System.reaction(2).product = [A0000,T2];
System.reaction(2).propensity = koff2 * A0020;

System.reaction(3).educt = A0200;
System.reaction(3).product = [A0000,T2];
System.reaction(3).propensity = koff2 * A0200;

System.reaction(4).educt = A2000;
System.reaction(4).product = [A0000,T2];
System.reaction(4).propensity = koff2 * A2000;

System.reaction(5).educt = [A0000,T2];
System.reaction(5).product = A0002;
System.reaction(5).propensity = kon2 * A0000 * T2;

System.reaction(6).educt = A0022;
System.reaction(6).product = [A0002,T2];
System.reaction(6).propensity = koff2 * A0022;

System.reaction(7).educt = A0202;
System.reaction(7).product = [A0002,T2];
System.reaction(7).propensity = koff2 * A0202;

System.reaction(8).educt = A2002;
System.reaction(8).product = [A0002,T2];
System.reaction(8).propensity = koff2 * A2002;

System.reaction(9).educt = A0002;
System.reaction(9).product = [A0002,R];
System.reaction(9).propensity = (r2) * A0002;

System.reaction(10).educt = [A0000,T2];
System.reaction(10).product = A0020;
System.reaction(10).propensity = kon2 * A0000 * T2;

System.reaction(11).educt = A0022;
System.reaction(11).product = [A0020,T2];
System.reaction(11).propensity = koff2 * A0022;

System.reaction(12).educt = A0220;
System.reaction(12).product = [A0020,T2];
System.reaction(12).propensity = koff2 * A0220;

System.reaction(13).educt = A2020;
System.reaction(13).product = [A0020,T2];
System.reaction(13).propensity = koff2 * A2020;

System.reaction(14).educt = A0020;
System.reaction(14).product = [A0020,R];
System.reaction(14).propensity = (r2) * A0020;

System.reaction(15).educt = [A0000,T2];
System.reaction(15).product = A0200;
System.reaction(15).propensity = kon2 * A0000 * T2;

System.reaction(16).educt = A0202;
System.reaction(16).product = [A0200,T2];
System.reaction(16).propensity = koff2 * A0202;

System.reaction(17).educt = A0220;
System.reaction(17).product = [A0200,T2];
System.reaction(17).propensity = koff2 * A0220;

System.reaction(18).educt = A2200;
System.reaction(18).product = [A0200,T2];
System.reaction(18).propensity = koff2 * A2200;

System.reaction(19).educt = A0200;
System.reaction(19).product = [A0200,R];
System.reaction(19).propensity = (r2) * A0200;

System.reaction(20).educt = [A0000,T2];
System.reaction(20).product = A2000;
System.reaction(20).propensity = kon2 * A0000 * T2;

System.reaction(21).educt = A2002;
System.reaction(21).product = [A2000,T2];
System.reaction(21).propensity = koff2 * A2002;

System.reaction(22).educt = A2020;
System.reaction(22).product = [A2000,T2];
System.reaction(22).propensity = koff2 * A2020;

System.reaction(23).educt = A2200;
System.reaction(23).product = [A2000,T2];
System.reaction(23).propensity = koff2 * A2200;

System.reaction(24).educt = A2000;
System.reaction(24).product = [A2000,R];
System.reaction(24).propensity = (r2) * A2000;

System.reaction(25).educt = [A0002,T2];
System.reaction(25).product = A0022;
System.reaction(25).propensity = kon2 * A0002 * T2;

System.reaction(26).educt = [A0020,T2];
System.reaction(26).product = A0022;
System.reaction(26).propensity = kon2 * A0020 * T2;

System.reaction(27).educt = A0222;
System.reaction(27).product = [A0022,T2];
System.reaction(27).propensity = koff2 * A0222;

System.reaction(28).educt = A2022;
System.reaction(28).product = [A0022,T2];
System.reaction(28).propensity = koff2 * A2022;

System.reaction(29).educt = A0022;
System.reaction(29).product = [A0022,R];
System.reaction(29).propensity = (r2 + r2) * A0022;

System.reaction(30).educt = [A0002,T2];
System.reaction(30).product = A0202;
System.reaction(30).propensity = kon2 * A0002 * T2;

System.reaction(31).educt = [A0200,T2];
System.reaction(31).product = A0202;
System.reaction(31).propensity = kon2 * A0200 * T2;

System.reaction(32).educt = A0222;
System.reaction(32).product = [A0202,T2];
System.reaction(32).propensity = koff2 * A0222;

System.reaction(33).educt = A2202;
System.reaction(33).product = [A0202,T2];
System.reaction(33).propensity = koff2 * A2202;

System.reaction(34).educt = A0202;
System.reaction(34).product = [A0202,R];
System.reaction(34).propensity = (r2 + r2) * A0202;

System.reaction(35).educt = [A0020,T2];
System.reaction(35).product = A0220;
System.reaction(35).propensity = kon2 * A0020 * T2;

System.reaction(36).educt = [A0200,T2];
System.reaction(36).product = A0220;
System.reaction(36).propensity = kon2 * A0200 * T2;

System.reaction(37).educt = A0222;
System.reaction(37).product = [A0220,T2];
System.reaction(37).propensity = koff2 * A0222;

System.reaction(38).educt = A2220;
System.reaction(38).product = [A0220,T2];
System.reaction(38).propensity = koff2 * A2220;

System.reaction(39).educt = A0220;
System.reaction(39).product = [A0220,R];
System.reaction(39).propensity = (r2 + r2) * A0220;

System.reaction(40).educt = [A0002,T2];
System.reaction(40).product = A2002;
System.reaction(40).propensity = kon2 * A0002 * T2;

System.reaction(41).educt = [A2000,T2];
System.reaction(41).product = A2002;
System.reaction(41).propensity = kon2 * A2000 * T2;

System.reaction(42).educt = A2022;
System.reaction(42).product = [A2002,T2];
System.reaction(42).propensity = koff2 * A2022;

System.reaction(43).educt = A2202;
System.reaction(43).product = [A2002,T2];
System.reaction(43).propensity = koff2 * A2202;

System.reaction(44).educt = A2002;
System.reaction(44).product = [A2002,R];
System.reaction(44).propensity = (r2 + r2) * A2002;

System.reaction(45).educt = [A0020,T2];
System.reaction(45).product = A2020;
System.reaction(45).propensity = kon2 * A0020 * T2;

System.reaction(46).educt = [A2000,T2];
System.reaction(46).product = A2020;
System.reaction(46).propensity = kon2 * A2000 * T2;

System.reaction(47).educt = A2022;
System.reaction(47).product = [A2020,T2];
System.reaction(47).propensity = koff2 * A2022;

System.reaction(48).educt = A2220;
System.reaction(48).product = [A2020,T2];
System.reaction(48).propensity = koff2 * A2220;

System.reaction(49).educt = A2020;
System.reaction(49).product = [A2020,R];
System.reaction(49).propensity = (r2 + r2) * A2020;

System.reaction(50).educt = [A0200,T2];
System.reaction(50).product = A2200;
System.reaction(50).propensity = kon2 * A0200 * T2;

System.reaction(51).educt = [A2000,T2];
System.reaction(51).product = A2200;
System.reaction(51).propensity = kon2 * A2000 * T2;

System.reaction(52).educt = A2202;
System.reaction(52).product = [A2200,T2];
System.reaction(52).propensity = koff2 * A2202;

System.reaction(53).educt = A2220;
System.reaction(53).product = [A2200,T2];
System.reaction(53).propensity = koff2 * A2220;

System.reaction(54).educt = A2200;
System.reaction(54).product = [A2200,R];
System.reaction(54).propensity = (r2 + r2) * A2200;

System.reaction(55).educt = [A0022,T2];
System.reaction(55).product = A0222;
System.reaction(55).propensity = kon2 * A0022 * T2;

System.reaction(56).educt = [A0202,T2];
System.reaction(56).product = A0222;
System.reaction(56).propensity = kon2 * A0202 * T2;

System.reaction(57).educt = [A0220,T2];
System.reaction(57).product = A0222;
System.reaction(57).propensity = kon2 * A0220 * T2;

System.reaction(58).educt = A2222;
System.reaction(58).product = [A0222,T2];
System.reaction(58).propensity = koff2 * A2222;

System.reaction(59).educt = A0222;
System.reaction(59).product = [A0222,R];
System.reaction(59).propensity = (r2 + r2 + r2) * A0222;

System.reaction(60).educt = [A0022,T2];
System.reaction(60).product = A2022;
System.reaction(60).propensity = kon2 * A0022 * T2;

System.reaction(61).educt = [A2002,T2];
System.reaction(61).product = A2022;
System.reaction(61).propensity = kon2 * A2002 * T2;

System.reaction(62).educt = [A2020,T2];
System.reaction(62).product = A2022;
System.reaction(62).propensity = kon2 * A2020 * T2;

System.reaction(63).educt = A2222;
System.reaction(63).product = [A2022,T2];
System.reaction(63).propensity = koff2 * A2222;

System.reaction(64).educt = A2022;
System.reaction(64).product = [A2022,R];
System.reaction(64).propensity = (r2 + r2 + r2) * A2022;

System.reaction(65).educt = [A0202,T2];
System.reaction(65).product = A2202;
System.reaction(65).propensity = kon2 * A0202 * T2;

System.reaction(66).educt = [A2002,T2];
System.reaction(66).product = A2202;
System.reaction(66).propensity = kon2 * A2002 * T2;

System.reaction(67).educt = [A2200,T2];
System.reaction(67).product = A2202;
System.reaction(67).propensity = kon2 * A2200 * T2;

System.reaction(68).educt = A2222;
System.reaction(68).product = [A2202,T2];
System.reaction(68).propensity = koff2 * A2222;

System.reaction(69).educt = A2202;
System.reaction(69).product = [A2202,R];
System.reaction(69).propensity = (r2 + r2 + r2) * A2202;

System.reaction(70).educt = [A0220,T2];
System.reaction(70).product = A2220;
System.reaction(70).propensity = kon2 * A0220 * T2;

System.reaction(71).educt = [A2020,T2];
System.reaction(71).product = A2220;
System.reaction(71).propensity = kon2 * A2020 * T2;

System.reaction(72).educt = [A2200,T2];
System.reaction(72).product = A2220;
System.reaction(72).propensity = kon2 * A2200 * T2;

System.reaction(73).educt = A2222;
System.reaction(73).product = [A2220,T2];
System.reaction(73).propensity = koff2 * A2222;

System.reaction(74).educt = A2220;
System.reaction(74).product = [A2220,R];
System.reaction(74).propensity = (r2 + r2 + r2) * A2220;

System.reaction(75).educt = [A0222,T2];
System.reaction(75).product = A2222;
System.reaction(75).propensity = kon2 * A0222 * T2;

System.reaction(76).educt = [A2022,T2];
System.reaction(76).product = A2222;
System.reaction(76).propensity = kon2 * A2022 * T2;

System.reaction(77).educt = [A2202,T2];
System.reaction(77).product = A2222;
System.reaction(77).propensity = kon2 * A2202 * T2;

System.reaction(78).educt = [A2220,T2];
System.reaction(78).product = A2222;
System.reaction(78).propensity = kon2 * A2220 * T2;

System.reaction(79).educt = A2222;
System.reaction(79).product = [A2222,R];
System.reaction(79).propensity = (r2 + r2 + r2 + r2) * A2222;

System.reaction(80).educt = T2;
System.reaction(80).product = [];
System.reaction(80).propensity = betam2 * T2;

System.reaction(81).educt = [];
System.reaction(81).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(81).propensity = beta2;

System.reaction(82).educt = R;
System.reaction(82).product = [];
System.reaction(82).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
