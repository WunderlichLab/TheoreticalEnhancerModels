syms A000 A002 A020 A022 A200 A202 A220 A222 T2 kon2 koff2 r2 myalpha beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A000; A002; A020; A022; A200; A202; A220; A222; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4) + x(5) + x(6) + x(7) + x(8)) == 1);
System.parameter.variable = [kon2; koff2; r2; myalpha; beta2; betam2];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A002;
System.reaction(1).product = [A000,T2];
System.reaction(1).propensity = koff2 * A002;

System.reaction(2).educt = A020;
System.reaction(2).product = [A000,T2];
System.reaction(2).propensity = koff2 * A020;

System.reaction(3).educt = A200;
System.reaction(3).product = [A000,T2];
System.reaction(3).propensity = koff2 * A200;

System.reaction(4).educt = [A000,T2];
System.reaction(4).product = A002;
System.reaction(4).propensity = kon2 * A000 * T2;

System.reaction(5).educt = A022;
System.reaction(5).product = [A002,T2];
System.reaction(5).propensity = koff2 * A022;

System.reaction(6).educt = A202;
System.reaction(6).product = [A002,T2];
System.reaction(6).propensity = koff2 * A202;

System.reaction(7).educt = A002;
System.reaction(7).product = [A002,R];
System.reaction(7).propensity = (r2) * A002;

System.reaction(8).educt = [A000,T2];
System.reaction(8).product = A020;
System.reaction(8).propensity = kon2 * A000 * T2;

System.reaction(9).educt = A022;
System.reaction(9).product = [A020,T2];
System.reaction(9).propensity = koff2 * A022;

System.reaction(10).educt = A220;
System.reaction(10).product = [A020,T2];
System.reaction(10).propensity = koff2 * A220;

System.reaction(11).educt = A020;
System.reaction(11).product = [A020,R];
System.reaction(11).propensity = (r2) * A020;

System.reaction(12).educt = [A000,T2];
System.reaction(12).product = A200;
System.reaction(12).propensity = kon2 * A000 * T2;

System.reaction(13).educt = A202;
System.reaction(13).product = [A200,T2];
System.reaction(13).propensity = koff2 * A202;

System.reaction(14).educt = A220;
System.reaction(14).product = [A200,T2];
System.reaction(14).propensity = koff2 * A220;

System.reaction(15).educt = A200;
System.reaction(15).product = [A200,R];
System.reaction(15).propensity = (r2) * A200;

System.reaction(16).educt = [A002,T2];
System.reaction(16).product = A022;
System.reaction(16).propensity = kon2 * A002 * T2;

System.reaction(17).educt = [A020,T2];
System.reaction(17).product = A022;
System.reaction(17).propensity = kon2 * A020 * T2;

System.reaction(18).educt = A222;
System.reaction(18).product = [A022,T2];
System.reaction(18).propensity = koff2 * A222;

System.reaction(19).educt = A022;
System.reaction(19).product = [A022,R];
System.reaction(19).propensity = (r2 + r2) * A022;

System.reaction(20).educt = [A002,T2];
System.reaction(20).product = A202;
System.reaction(20).propensity = kon2 * A002 * T2;

System.reaction(21).educt = [A200,T2];
System.reaction(21).product = A202;
System.reaction(21).propensity = kon2 * A200 * T2;

System.reaction(22).educt = A222;
System.reaction(22).product = [A202,T2];
System.reaction(22).propensity = koff2 * A222;

System.reaction(23).educt = A202;
System.reaction(23).product = [A202,R];
System.reaction(23).propensity = (r2 + r2) * A202;

System.reaction(24).educt = [A020,T2];
System.reaction(24).product = A220;
System.reaction(24).propensity = kon2 * A020 * T2;

System.reaction(25).educt = [A200,T2];
System.reaction(25).product = A220;
System.reaction(25).propensity = kon2 * A200 * T2;

System.reaction(26).educt = A222;
System.reaction(26).product = [A220,T2];
System.reaction(26).propensity = koff2 * A222;

System.reaction(27).educt = A220;
System.reaction(27).product = [A220,R];
System.reaction(27).propensity = (r2 + r2) * A220;

System.reaction(28).educt = [A022,T2];
System.reaction(28).product = A222;
System.reaction(28).propensity = kon2 * A022 * T2;

System.reaction(29).educt = [A202,T2];
System.reaction(29).product = A222;
System.reaction(29).propensity = kon2 * A202 * T2;

System.reaction(30).educt = [A220,T2];
System.reaction(30).product = A222;
System.reaction(30).propensity = kon2 * A220 * T2;

System.reaction(31).educt = A222;
System.reaction(31).product = [A222,R];
System.reaction(31).propensity = (r2 + r2 + r2) * A222;

System.reaction(32).educt = T2;
System.reaction(32).product = [];
System.reaction(32).propensity = betam2 * T2;

System.reaction(33).educt = [];
System.reaction(33).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(33).propensity = beta2;

System.reaction(34).educt = R;
System.reaction(34).product = [];
System.reaction(34).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
