syms A0 A2 B0 B2 T2 kon2 koff2 r2 myalpha beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0; A2; B0; B2; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 1; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2)) == 1 && (x(3) + x(4)) == 1);
System.parameter.variable = [kon2; koff2; r2; myalpha; beta2; betam2];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A2;
System.reaction(1).product = [A0,T2];
System.reaction(1).propensity = koff2 * A2;

System.reaction(2).educt = [A0,T2];
System.reaction(2).product = A2;
System.reaction(2).propensity = kon2 * A0 * T2;

System.reaction(3).educt = A2;
System.reaction(3).product = [A2,R];
System.reaction(3).propensity = (r2) * A2;

System.reaction(4).educt = B2;
System.reaction(4).product = [B0,T2];
System.reaction(4).propensity = koff2 * B2;

System.reaction(5).educt = [B0,T2];
System.reaction(5).product = B2;
System.reaction(5).propensity = kon2 * B0 * T2;

System.reaction(6).educt = B2;
System.reaction(6).product = [B2,R];
System.reaction(6).propensity = (r2) * B2;

System.reaction(7).educt = T2;
System.reaction(7).product = [];
System.reaction(7).propensity = betam2 * T2;

System.reaction(8).educt = [];
System.reaction(8).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(8).propensity = beta2;

System.reaction(9).educt = R;
System.reaction(9).product = [];
System.reaction(9).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
