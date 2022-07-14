syms A000 A001 A010 A011 A100 A101 A110 A111 T1 kon1 koff1 r1 myalpha beta1 betam1 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A000; A001; A010; A011; A100; A101; A110; A111; T1; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8)) == 1);
System.parameter.variable = [kon1; koff1; r1; myalpha; beta1; betam1];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A001;
System.reaction(1).product = [A000,T1];
System.reaction(1).propensity = koff1 * A001;

System.reaction(2).educt = A010;
System.reaction(2).product = [A000,T1];
System.reaction(2).propensity = koff1 * A010;

System.reaction(3).educt = A100;
System.reaction(3).product = [A000,T1];
System.reaction(3).propensity = koff1 * A100;

System.reaction(4).educt = [A000,T1];
System.reaction(4).product = A001;
System.reaction(4).propensity = kon1 * A000 * T1;

System.reaction(5).educt = A011;
System.reaction(5).product = [A001,T1];
System.reaction(5).propensity = koff1 * A011;

System.reaction(6).educt = A101;
System.reaction(6).product = [A001,T1];
System.reaction(6).propensity = koff1 * A101;

System.reaction(7).educt = A001;
System.reaction(7).product = [A001,R];
System.reaction(7).propensity = (r1) * A001;

System.reaction(8).educt = [A000,T1];
System.reaction(8).product = A010;
System.reaction(8).propensity = kon1 * A000 * T1;

System.reaction(9).educt = A011;
System.reaction(9).product = [A010,T1];
System.reaction(9).propensity = koff1 * A011;

System.reaction(10).educt = A110;
System.reaction(10).product = [A010,T1];
System.reaction(10).propensity = koff1 * A110;

System.reaction(11).educt = A010;
System.reaction(11).product = [A010,R];
System.reaction(11).propensity = (r1) * A010;

System.reaction(12).educt = [A000,T1];
System.reaction(12).product = A100;
System.reaction(12).propensity = kon1 * A000 * T1;

System.reaction(13).educt = A101;
System.reaction(13).product = [A100,T1];
System.reaction(13).propensity = koff1 * A101;

System.reaction(14).educt = A110;
System.reaction(14).product = [A100,T1];
System.reaction(14).propensity = koff1 * A110;

System.reaction(15).educt = A100;
System.reaction(15).product = [A100,R];
System.reaction(15).propensity = (r1) * A100;

System.reaction(16).educt = [A001,T1];
System.reaction(16).product = A011;
System.reaction(16).propensity = kon1 * A001 * T1;

System.reaction(17).educt = [A010,T1];
System.reaction(17).product = A011;
System.reaction(17).propensity = kon1 * A010 * T1;

System.reaction(18).educt = A111;
System.reaction(18).product = [A011,T1];
System.reaction(18).propensity = koff1 * A111;

System.reaction(19).educt = A011;
System.reaction(19).product = [A011,R];
System.reaction(19).propensity = (r1 + r1) * A011;

System.reaction(20).educt = [A001,T1];
System.reaction(20).product = A101;
System.reaction(20).propensity = kon1 * A001 * T1;

System.reaction(21).educt = [A100,T1];
System.reaction(21).product = A101;
System.reaction(21).propensity = kon1 * A100 * T1;

System.reaction(22).educt = A111;
System.reaction(22).product = [A101,T1];
System.reaction(22).propensity = koff1 * A111;

System.reaction(23).educt = A101;
System.reaction(23).product = [A101,R];
System.reaction(23).propensity = (r1 + r1) * A101;

System.reaction(24).educt = [A010,T1];
System.reaction(24).product = A110;
System.reaction(24).propensity = kon1 * A010 * T1;

System.reaction(25).educt = [A100,T1];
System.reaction(25).product = A110;
System.reaction(25).propensity = kon1 * A100 * T1;

System.reaction(26).educt = A111;
System.reaction(26).product = [A110,T1];
System.reaction(26).propensity = koff1 * A111;

System.reaction(27).educt = A110;
System.reaction(27).product = [A110,R];
System.reaction(27).propensity = (r1 + r1) * A110;

System.reaction(28).educt = [A011,T1];
System.reaction(28).product = A111;
System.reaction(28).propensity = kon1 * A011 * T1;

System.reaction(29).educt = [A101,T1];
System.reaction(29).product = A111;
System.reaction(29).propensity = kon1 * A101 * T1;

System.reaction(30).educt = [A110,T1];
System.reaction(30).product = A111;
System.reaction(30).propensity = kon1 * A110 * T1;

System.reaction(31).educt = A111;
System.reaction(31).product = [A111,R];
System.reaction(31).propensity = (r1 + r1 + r1) * A111;

System.reaction(32).educt = T1;
System.reaction(32).product = [];
System.reaction(32).propensity = betam1 * T1;

System.reaction(33).educt = [];
System.reaction(33).product = [T1, T1, T1, T1];
System.reaction(33).propensity = beta1;

System.reaction(34).educt = R;
System.reaction(34).product = [];
System.reaction(34).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
