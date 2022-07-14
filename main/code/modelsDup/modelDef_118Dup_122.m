syms A000 A002 A020 A022 A100 A102 A120 A122 B000 B002 B020 B022 B100 B102 B120 B122 T1 T2 kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A000; A002; A020; A022; A100; A102; A120; A122; B000; B002; B020; B022; B100; B102; B120; B122; T1; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 10; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8)) == 1 && (x(9) + x(10) + x(11) + x(12) + x(13) + x(14) + x(15) + x(16)) == 1);
System.parameter.variable = [kon1; koff1; r1; kon2; koff2; r2; myalpha; beta1; betam1; beta2; betam2];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A002;
System.reaction(1).product = [A000,T2];
System.reaction(1).propensity = koff2 * A002;

System.reaction(2).educt = A020;
System.reaction(2).product = [A000,T2];
System.reaction(2).propensity = koff2 * A020;

System.reaction(3).educt = A100;
System.reaction(3).product = [A000,T1];
System.reaction(3).propensity = koff1 * A100;

System.reaction(4).educt = [A000,T2];
System.reaction(4).product = A002;
System.reaction(4).propensity = kon2 * A000 * T2;

System.reaction(5).educt = A022;
System.reaction(5).product = [A002,T2];
System.reaction(5).propensity = koff2 * A022;

System.reaction(6).educt = A102;
System.reaction(6).product = [A002,T1];
System.reaction(6).propensity = koff1 * A102;

System.reaction(7).educt = A002;
System.reaction(7).product = [A002,R];
System.reaction(7).propensity = (r2) * A002;

System.reaction(8).educt = [A000,T2];
System.reaction(8).product = A020;
System.reaction(8).propensity = kon2 * A000 * T2;

System.reaction(9).educt = A022;
System.reaction(9).product = [A020,T2];
System.reaction(9).propensity = koff2 * A022;

System.reaction(10).educt = A120;
System.reaction(10).product = [A020,T1];
System.reaction(10).propensity = koff1 * A120;

System.reaction(11).educt = A020;
System.reaction(11).product = [A020,R];
System.reaction(11).propensity = (r2) * A020;

System.reaction(12).educt = [A000,T1];
System.reaction(12).product = A100;
System.reaction(12).propensity = kon1 * A000 * T1;

System.reaction(13).educt = A102;
System.reaction(13).product = [A100,T2];
System.reaction(13).propensity = koff2 * A102;

System.reaction(14).educt = A120;
System.reaction(14).product = [A100,T2];
System.reaction(14).propensity = koff2 * A120;

System.reaction(15).educt = A100;
System.reaction(15).product = [A100,R];
System.reaction(15).propensity = (r1) * A100;

System.reaction(16).educt = [A002,T2];
System.reaction(16).product = A022;
System.reaction(16).propensity = kon2 * A002 * T2;

System.reaction(17).educt = [A020,T2];
System.reaction(17).product = A022;
System.reaction(17).propensity = kon2 * A020 * T2;

System.reaction(18).educt = A122;
System.reaction(18).product = [A022,T1];
System.reaction(18).propensity = koff1 * A122;

System.reaction(19).educt = A022;
System.reaction(19).product = [A022,R];
System.reaction(19).propensity = (r2 + r2) * A022;

System.reaction(20).educt = [A002,T1];
System.reaction(20).product = A102;
System.reaction(20).propensity = kon1 * A002 * T1;

System.reaction(21).educt = [A100,T2];
System.reaction(21).product = A102;
System.reaction(21).propensity = kon2 * A100 * T2;

System.reaction(22).educt = A122;
System.reaction(22).product = [A102,T2];
System.reaction(22).propensity = koff2 * A122;

System.reaction(23).educt = A102;
System.reaction(23).product = [A102,R];
System.reaction(23).propensity = (r1 + r2) * A102;

System.reaction(24).educt = [A020,T1];
System.reaction(24).product = A120;
System.reaction(24).propensity = kon1 * A020 * T1;

System.reaction(25).educt = [A100,T2];
System.reaction(25).product = A120;
System.reaction(25).propensity = kon2 * A100 * T2;

System.reaction(26).educt = A122;
System.reaction(26).product = [A120,T2];
System.reaction(26).propensity = koff2 * A122;

System.reaction(27).educt = A120;
System.reaction(27).product = [A120,R];
System.reaction(27).propensity = (r1 + r2) * A120;

System.reaction(28).educt = [A022,T1];
System.reaction(28).product = A122;
System.reaction(28).propensity = kon1 * A022 * T1;

System.reaction(29).educt = [A102,T2];
System.reaction(29).product = A122;
System.reaction(29).propensity = kon2 * A102 * T2;

System.reaction(30).educt = [A120,T2];
System.reaction(30).product = A122;
System.reaction(30).propensity = kon2 * A120 * T2;

System.reaction(31).educt = A122;
System.reaction(31).product = [A122,R];
System.reaction(31).propensity = (r1 + r2 + r2) * A122;

System.reaction(32).educt = B002;
System.reaction(32).product = [B000,T2];
System.reaction(32).propensity = koff2 * B002;

System.reaction(33).educt = B020;
System.reaction(33).product = [B000,T2];
System.reaction(33).propensity = koff2 * B020;

System.reaction(34).educt = B100;
System.reaction(34).product = [B000,T1];
System.reaction(34).propensity = koff1 * B100;

System.reaction(35).educt = [B000,T2];
System.reaction(35).product = B002;
System.reaction(35).propensity = kon2 * B000 * T2;

System.reaction(36).educt = B022;
System.reaction(36).product = [B002,T2];
System.reaction(36).propensity = koff2 * B022;

System.reaction(37).educt = B102;
System.reaction(37).product = [B002,T1];
System.reaction(37).propensity = koff1 * B102;

System.reaction(38).educt = B002;
System.reaction(38).product = [B002,R];
System.reaction(38).propensity = (r2) * B002;

System.reaction(39).educt = [B000,T2];
System.reaction(39).product = B020;
System.reaction(39).propensity = kon2 * B000 * T2;

System.reaction(40).educt = B022;
System.reaction(40).product = [B020,T2];
System.reaction(40).propensity = koff2 * B022;

System.reaction(41).educt = B120;
System.reaction(41).product = [B020,T1];
System.reaction(41).propensity = koff1 * B120;

System.reaction(42).educt = B020;
System.reaction(42).product = [B020,R];
System.reaction(42).propensity = (r2) * B020;

System.reaction(43).educt = [B000,T1];
System.reaction(43).product = B100;
System.reaction(43).propensity = kon1 * B000 * T1;

System.reaction(44).educt = B102;
System.reaction(44).product = [B100,T2];
System.reaction(44).propensity = koff2 * B102;

System.reaction(45).educt = B120;
System.reaction(45).product = [B100,T2];
System.reaction(45).propensity = koff2 * B120;

System.reaction(46).educt = B100;
System.reaction(46).product = [B100,R];
System.reaction(46).propensity = (r1) * B100;

System.reaction(47).educt = [B002,T2];
System.reaction(47).product = B022;
System.reaction(47).propensity = kon2 * B002 * T2;

System.reaction(48).educt = [B020,T2];
System.reaction(48).product = B022;
System.reaction(48).propensity = kon2 * B020 * T2;

System.reaction(49).educt = B122;
System.reaction(49).product = [B022,T1];
System.reaction(49).propensity = koff1 * B122;

System.reaction(50).educt = B022;
System.reaction(50).product = [B022,R];
System.reaction(50).propensity = (r2 + r2) * B022;

System.reaction(51).educt = [B002,T1];
System.reaction(51).product = B102;
System.reaction(51).propensity = kon1 * B002 * T1;

System.reaction(52).educt = [B100,T2];
System.reaction(52).product = B102;
System.reaction(52).propensity = kon2 * B100 * T2;

System.reaction(53).educt = B122;
System.reaction(53).product = [B102,T2];
System.reaction(53).propensity = koff2 * B122;

System.reaction(54).educt = B102;
System.reaction(54).product = [B102,R];
System.reaction(54).propensity = (r1 + r2) * B102;

System.reaction(55).educt = [B020,T1];
System.reaction(55).product = B120;
System.reaction(55).propensity = kon1 * B020 * T1;

System.reaction(56).educt = [B100,T2];
System.reaction(56).product = B120;
System.reaction(56).propensity = kon2 * B100 * T2;

System.reaction(57).educt = B122;
System.reaction(57).product = [B120,T2];
System.reaction(57).propensity = koff2 * B122;

System.reaction(58).educt = B120;
System.reaction(58).product = [B120,R];
System.reaction(58).propensity = (r1 + r2) * B120;

System.reaction(59).educt = [B022,T1];
System.reaction(59).product = B122;
System.reaction(59).propensity = kon1 * B022 * T1;

System.reaction(60).educt = [B102,T2];
System.reaction(60).product = B122;
System.reaction(60).propensity = kon2 * B102 * T2;

System.reaction(61).educt = [B120,T2];
System.reaction(61).product = B122;
System.reaction(61).propensity = kon2 * B120 * T2;

System.reaction(62).educt = B122;
System.reaction(62).product = [B122,R];
System.reaction(62).propensity = (r1 + r2 + r2) * B122;

System.reaction(63).educt = T1;
System.reaction(63).product = [];
System.reaction(63).propensity = betam1 * T1;

System.reaction(64).educt = [];
System.reaction(64).product = [T1, T1, T1, T1];
System.reaction(64).propensity = beta1;

System.reaction(65).educt = T2;
System.reaction(65).product = [];
System.reaction(65).propensity = betam2 * T2;

System.reaction(66).educt = [];
System.reaction(66).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(66).propensity = beta2;

System.reaction(67).educt = R;
System.reaction(67).product = [];
System.reaction(67).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
