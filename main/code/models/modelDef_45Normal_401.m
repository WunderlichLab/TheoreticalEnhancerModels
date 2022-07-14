syms A0000 A0001 A0010 A0011 A0100 A0101 A0110 A0111 A1000 A1001 A1010 A1011 A1100 A1101 A1110 A1111 T1 kon1 koff1 r1 myalpha beta1 betam1 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0000; A0001; A0010; A0011; A0100; A0101; A0110; A0111; A1000; A1001; A1010; A1011; A1100; A1101; A1110; A1111; T1; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8) + x(9) + x(10) + x(11) + x(12) + x(13) + x(14) + x(15) + x(16)) == 1);
System.parameter.variable = [kon1; koff1; r1; myalpha; beta1; betam1];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A0001;
System.reaction(1).product = [A0000,T1];
System.reaction(1).propensity = koff1 * A0001;

System.reaction(2).educt = A0010;
System.reaction(2).product = [A0000,T1];
System.reaction(2).propensity = koff1 * A0010;

System.reaction(3).educt = A0100;
System.reaction(3).product = [A0000,T1];
System.reaction(3).propensity = koff1 * A0100;

System.reaction(4).educt = A1000;
System.reaction(4).product = [A0000,T1];
System.reaction(4).propensity = koff1 * A1000;

System.reaction(5).educt = [A0000,T1];
System.reaction(5).product = A0001;
System.reaction(5).propensity = kon1 * A0000 * T1;

System.reaction(6).educt = A0011;
System.reaction(6).product = [A0001,T1];
System.reaction(6).propensity = koff1 * A0011;

System.reaction(7).educt = A0101;
System.reaction(7).product = [A0001,T1];
System.reaction(7).propensity = koff1 * A0101;

System.reaction(8).educt = A1001;
System.reaction(8).product = [A0001,T1];
System.reaction(8).propensity = koff1 * A1001;

System.reaction(9).educt = A0001;
System.reaction(9).product = [A0001,R];
System.reaction(9).propensity = (r1) * A0001;

System.reaction(10).educt = [A0000,T1];
System.reaction(10).product = A0010;
System.reaction(10).propensity = kon1 * A0000 * T1;

System.reaction(11).educt = A0011;
System.reaction(11).product = [A0010,T1];
System.reaction(11).propensity = koff1 * A0011;

System.reaction(12).educt = A0110;
System.reaction(12).product = [A0010,T1];
System.reaction(12).propensity = koff1 * A0110;

System.reaction(13).educt = A1010;
System.reaction(13).product = [A0010,T1];
System.reaction(13).propensity = koff1 * A1010;

System.reaction(14).educt = A0010;
System.reaction(14).product = [A0010,R];
System.reaction(14).propensity = (r1) * A0010;

System.reaction(15).educt = [A0000,T1];
System.reaction(15).product = A0100;
System.reaction(15).propensity = kon1 * A0000 * T1;

System.reaction(16).educt = A0101;
System.reaction(16).product = [A0100,T1];
System.reaction(16).propensity = koff1 * A0101;

System.reaction(17).educt = A0110;
System.reaction(17).product = [A0100,T1];
System.reaction(17).propensity = koff1 * A0110;

System.reaction(18).educt = A1100;
System.reaction(18).product = [A0100,T1];
System.reaction(18).propensity = koff1 * A1100;

System.reaction(19).educt = A0100;
System.reaction(19).product = [A0100,R];
System.reaction(19).propensity = (r1) * A0100;

System.reaction(20).educt = [A0000,T1];
System.reaction(20).product = A1000;
System.reaction(20).propensity = kon1 * A0000 * T1;

System.reaction(21).educt = A1001;
System.reaction(21).product = [A1000,T1];
System.reaction(21).propensity = koff1 * A1001;

System.reaction(22).educt = A1010;
System.reaction(22).product = [A1000,T1];
System.reaction(22).propensity = koff1 * A1010;

System.reaction(23).educt = A1100;
System.reaction(23).product = [A1000,T1];
System.reaction(23).propensity = koff1 * A1100;

System.reaction(24).educt = A1000;
System.reaction(24).product = [A1000,R];
System.reaction(24).propensity = (r1) * A1000;

System.reaction(25).educt = [A0001,T1];
System.reaction(25).product = A0011;
System.reaction(25).propensity = kon1 * A0001 * T1;

System.reaction(26).educt = [A0010,T1];
System.reaction(26).product = A0011;
System.reaction(26).propensity = kon1 * A0010 * T1;

System.reaction(27).educt = A0111;
System.reaction(27).product = [A0011,T1];
System.reaction(27).propensity = koff1 * A0111;

System.reaction(28).educt = A1011;
System.reaction(28).product = [A0011,T1];
System.reaction(28).propensity = koff1 * A1011;

System.reaction(29).educt = A0011;
System.reaction(29).product = [A0011,R];
System.reaction(29).propensity = (r1 + r1) * A0011;

System.reaction(30).educt = [A0001,T1];
System.reaction(30).product = A0101;
System.reaction(30).propensity = kon1 * A0001 * T1;

System.reaction(31).educt = [A0100,T1];
System.reaction(31).product = A0101;
System.reaction(31).propensity = kon1 * A0100 * T1;

System.reaction(32).educt = A0111;
System.reaction(32).product = [A0101,T1];
System.reaction(32).propensity = koff1 * A0111;

System.reaction(33).educt = A1101;
System.reaction(33).product = [A0101,T1];
System.reaction(33).propensity = koff1 * A1101;

System.reaction(34).educt = A0101;
System.reaction(34).product = [A0101,R];
System.reaction(34).propensity = (r1 + r1) * A0101;

System.reaction(35).educt = [A0010,T1];
System.reaction(35).product = A0110;
System.reaction(35).propensity = kon1 * A0010 * T1;

System.reaction(36).educt = [A0100,T1];
System.reaction(36).product = A0110;
System.reaction(36).propensity = kon1 * A0100 * T1;

System.reaction(37).educt = A0111;
System.reaction(37).product = [A0110,T1];
System.reaction(37).propensity = koff1 * A0111;

System.reaction(38).educt = A1110;
System.reaction(38).product = [A0110,T1];
System.reaction(38).propensity = koff1 * A1110;

System.reaction(39).educt = A0110;
System.reaction(39).product = [A0110,R];
System.reaction(39).propensity = (r1 + r1) * A0110;

System.reaction(40).educt = [A0001,T1];
System.reaction(40).product = A1001;
System.reaction(40).propensity = kon1 * A0001 * T1;

System.reaction(41).educt = [A1000,T1];
System.reaction(41).product = A1001;
System.reaction(41).propensity = kon1 * A1000 * T1;

System.reaction(42).educt = A1011;
System.reaction(42).product = [A1001,T1];
System.reaction(42).propensity = koff1 * A1011;

System.reaction(43).educt = A1101;
System.reaction(43).product = [A1001,T1];
System.reaction(43).propensity = koff1 * A1101;

System.reaction(44).educt = A1001;
System.reaction(44).product = [A1001,R];
System.reaction(44).propensity = (r1 + r1) * A1001;

System.reaction(45).educt = [A0010,T1];
System.reaction(45).product = A1010;
System.reaction(45).propensity = kon1 * A0010 * T1;

System.reaction(46).educt = [A1000,T1];
System.reaction(46).product = A1010;
System.reaction(46).propensity = kon1 * A1000 * T1;

System.reaction(47).educt = A1011;
System.reaction(47).product = [A1010,T1];
System.reaction(47).propensity = koff1 * A1011;

System.reaction(48).educt = A1110;
System.reaction(48).product = [A1010,T1];
System.reaction(48).propensity = koff1 * A1110;

System.reaction(49).educt = A1010;
System.reaction(49).product = [A1010,R];
System.reaction(49).propensity = (r1 + r1) * A1010;

System.reaction(50).educt = [A0100,T1];
System.reaction(50).product = A1100;
System.reaction(50).propensity = kon1 * A0100 * T1;

System.reaction(51).educt = [A1000,T1];
System.reaction(51).product = A1100;
System.reaction(51).propensity = kon1 * A1000 * T1;

System.reaction(52).educt = A1101;
System.reaction(52).product = [A1100,T1];
System.reaction(52).propensity = koff1 * A1101;

System.reaction(53).educt = A1110;
System.reaction(53).product = [A1100,T1];
System.reaction(53).propensity = koff1 * A1110;

System.reaction(54).educt = A1100;
System.reaction(54).product = [A1100,R];
System.reaction(54).propensity = (r1 + r1) * A1100;

System.reaction(55).educt = [A0011,T1];
System.reaction(55).product = A0111;
System.reaction(55).propensity = kon1 * A0011 * T1;

System.reaction(56).educt = [A0101,T1];
System.reaction(56).product = A0111;
System.reaction(56).propensity = kon1 * A0101 * T1;

System.reaction(57).educt = [A0110,T1];
System.reaction(57).product = A0111;
System.reaction(57).propensity = kon1 * A0110 * T1;

System.reaction(58).educt = A1111;
System.reaction(58).product = [A0111,T1];
System.reaction(58).propensity = koff1 * A1111;

System.reaction(59).educt = A0111;
System.reaction(59).product = [A0111,R];
System.reaction(59).propensity = (r1 + r1 + r1) * A0111;

System.reaction(60).educt = [A0011,T1];
System.reaction(60).product = A1011;
System.reaction(60).propensity = kon1 * A0011 * T1;

System.reaction(61).educt = [A1001,T1];
System.reaction(61).product = A1011;
System.reaction(61).propensity = kon1 * A1001 * T1;

System.reaction(62).educt = [A1010,T1];
System.reaction(62).product = A1011;
System.reaction(62).propensity = kon1 * A1010 * T1;

System.reaction(63).educt = A1111;
System.reaction(63).product = [A1011,T1];
System.reaction(63).propensity = koff1 * A1111;

System.reaction(64).educt = A1011;
System.reaction(64).product = [A1011,R];
System.reaction(64).propensity = (r1 + r1 + r1) * A1011;

System.reaction(65).educt = [A0101,T1];
System.reaction(65).product = A1101;
System.reaction(65).propensity = kon1 * A0101 * T1;

System.reaction(66).educt = [A1001,T1];
System.reaction(66).product = A1101;
System.reaction(66).propensity = kon1 * A1001 * T1;

System.reaction(67).educt = [A1100,T1];
System.reaction(67).product = A1101;
System.reaction(67).propensity = kon1 * A1100 * T1;

System.reaction(68).educt = A1111;
System.reaction(68).product = [A1101,T1];
System.reaction(68).propensity = koff1 * A1111;

System.reaction(69).educt = A1101;
System.reaction(69).product = [A1101,R];
System.reaction(69).propensity = (r1 + r1 + r1) * A1101;

System.reaction(70).educt = [A0110,T1];
System.reaction(70).product = A1110;
System.reaction(70).propensity = kon1 * A0110 * T1;

System.reaction(71).educt = [A1010,T1];
System.reaction(71).product = A1110;
System.reaction(71).propensity = kon1 * A1010 * T1;

System.reaction(72).educt = [A1100,T1];
System.reaction(72).product = A1110;
System.reaction(72).propensity = kon1 * A1100 * T1;

System.reaction(73).educt = A1111;
System.reaction(73).product = [A1110,T1];
System.reaction(73).propensity = koff1 * A1111;

System.reaction(74).educt = A1110;
System.reaction(74).product = [A1110,R];
System.reaction(74).propensity = (r1 + r1 + r1) * A1110;

System.reaction(75).educt = [A0111,T1];
System.reaction(75).product = A1111;
System.reaction(75).propensity = kon1 * A0111 * T1;

System.reaction(76).educt = [A1011,T1];
System.reaction(76).product = A1111;
System.reaction(76).propensity = kon1 * A1011 * T1;

System.reaction(77).educt = [A1101,T1];
System.reaction(77).product = A1111;
System.reaction(77).propensity = kon1 * A1101 * T1;

System.reaction(78).educt = [A1110,T1];
System.reaction(78).product = A1111;
System.reaction(78).propensity = kon1 * A1110 * T1;

System.reaction(79).educt = A1111;
System.reaction(79).product = [A1111,R];
System.reaction(79).propensity = (r1 + r1 + r1 + r1) * A1111;

System.reaction(80).educt = T1;
System.reaction(80).product = [];
System.reaction(80).propensity = betam1 * T1;

System.reaction(81).educt = [];
System.reaction(81).product = [T1, T1, T1, T1];
System.reaction(81).propensity = beta1;

System.reaction(82).educt = R;
System.reaction(82).product = [];
System.reaction(82).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
