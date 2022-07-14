syms A0 A1 B0 B1 C0 C2 T1 T2 kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0; A1; B0; B1; C0; C2; T1; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 10; 10; 150];
System.state.mu0 = [1; 0; 1; 0; 1; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2)) == 1 && (x(3) + x(4)) == 1 && (x(5) + x(6)) == 1);
System.parameter.variable = [kon1; koff1; r1; kon2; koff2; r2; myalpha; beta1; betam1; beta2; betam2];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A1;
System.reaction(1).product = [A0,T1];
System.reaction(1).propensity = koff1 * A1;

System.reaction(2).educt = [A0,T1];
System.reaction(2).product = A1;
System.reaction(2).propensity = kon1 * A0 * T1;

System.reaction(3).educt = A1;
System.reaction(3).product = [A1,R];
System.reaction(3).propensity = (r1) * A1;

System.reaction(4).educt = B1;
System.reaction(4).product = [B0,T1];
System.reaction(4).propensity = koff1 * B1;

System.reaction(5).educt = [B0,T1];
System.reaction(5).product = B1;
System.reaction(5).propensity = kon1 * B0 * T1;

System.reaction(6).educt = B1;
System.reaction(6).product = [B1,R];
System.reaction(6).propensity = (r1) * B1;

System.reaction(7).educt = C2;
System.reaction(7).product = [C0,T2];
System.reaction(7).propensity = koff2 * C2;

System.reaction(8).educt = [C0,T2];
System.reaction(8).product = C2;
System.reaction(8).propensity = kon2 * C0 * T2;

System.reaction(9).educt = C2;
System.reaction(9).product = [C2,R];
System.reaction(9).propensity = (r2) * C2;

System.reaction(10).educt = T1;
System.reaction(10).product = [];
System.reaction(10).propensity = betam1 * T1;

System.reaction(11).educt = [];
System.reaction(11).product = [T1, T1, T1, T1];
System.reaction(11).propensity = beta1;

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
