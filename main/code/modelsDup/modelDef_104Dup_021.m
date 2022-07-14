syms A00 A02 A20 A22 T2 kon2 koff2 r2 myalpha beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A00; A02; A20; A22; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4)) == 1);
System.parameter.variable = [kon2; koff2; r2; myalpha; beta2; betam2];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A02;
System.reaction(1).product = [A00,T2];
System.reaction(1).propensity = koff2 * A02;

System.reaction(2).educt = A20;
System.reaction(2).product = [A00,T2];
System.reaction(2).propensity = koff2 * A20;

System.reaction(3).educt = [A00,T2];
System.reaction(3).product = A02;
System.reaction(3).propensity = kon2 * A00 * T2;

System.reaction(4).educt = A22;
System.reaction(4).product = [A02,T2];
System.reaction(4).propensity = koff2 * A22;

System.reaction(5).educt = A02;
System.reaction(5).product = [A02,R];
System.reaction(5).propensity = (r2) * A02;

System.reaction(6).educt = [A00,T2];
System.reaction(6).product = A20;
System.reaction(6).propensity = kon2 * A00 * T2;

System.reaction(7).educt = A22;
System.reaction(7).product = [A20,T2];
System.reaction(7).propensity = koff2 * A22;

System.reaction(8).educt = A20;
System.reaction(8).product = [A20,R];
System.reaction(8).propensity = (r2) * A20;

System.reaction(9).educt = [A02,T2];
System.reaction(9).product = A22;
System.reaction(9).propensity = kon2 * A02 * T2;

System.reaction(10).educt = [A20,T2];
System.reaction(10).product = A22;
System.reaction(10).propensity = kon2 * A20 * T2;

System.reaction(11).educt = A22;
System.reaction(11).product = [A22,R];
System.reaction(11).propensity = (r2 + r2) * A22;

System.reaction(12).educt = T2;
System.reaction(12).product = [];
System.reaction(12).propensity = betam2 * T2;

System.reaction(13).educt = [];
System.reaction(13).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(13).propensity = beta2;

System.reaction(14).educt = R;
System.reaction(14).product = [];
System.reaction(14).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
