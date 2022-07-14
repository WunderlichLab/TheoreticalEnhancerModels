syms A0 A1 T1 kon1 koff1 r1 myalpha beta1 betam1 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0; A1; T1; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0];
System.state.xmax = [1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2)) == 1);
System.parameter.variable = [kon1; koff1; r1; myalpha; beta1; betam1];
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

System.reaction(4).educt = T1;
System.reaction(4).product = [];
System.reaction(4).propensity = betam1 * T1;

System.reaction(5).educt = [];
System.reaction(5).product = [T1, T1, T1, T1];
System.reaction(5).propensity = beta1;

System.reaction(6).educt = R;
System.reaction(6).product = [];
System.reaction(6).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
