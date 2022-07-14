syms A00 A01 A10 A11 T1 kon1 koff1 r1 myalpha beta1 betam1 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A00; A01; A10; A11; T1; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4)) == 1);
System.parameter.variable = [kon1; koff1; r1; myalpha; beta1; betam1];
System.kappa.variable = [Omega];
System.scaleIndicator = 'microscopic';

System.reaction(1).educt = A01;
System.reaction(1).product = [A00,T1];
System.reaction(1).propensity = koff1 * A01;

System.reaction(2).educt = A10;
System.reaction(2).product = [A00,T1];
System.reaction(2).propensity = koff1 * A10;

System.reaction(3).educt = [A00,T1];
System.reaction(3).product = A01;
System.reaction(3).propensity = kon1 * A00 * T1;

System.reaction(4).educt = A11;
System.reaction(4).product = [A01,T1];
System.reaction(4).propensity = koff1 * A11;

System.reaction(5).educt = A01;
System.reaction(5).product = [A01,R];
System.reaction(5).propensity = (r1) * A01;

System.reaction(6).educt = [A00,T1];
System.reaction(6).product = A10;
System.reaction(6).propensity = kon1 * A00 * T1;

System.reaction(7).educt = A11;
System.reaction(7).product = [A10,T1];
System.reaction(7).propensity = koff1 * A11;

System.reaction(8).educt = A10;
System.reaction(8).product = [A10,R];
System.reaction(8).propensity = (r1) * A10;

System.reaction(9).educt = [A01,T1];
System.reaction(9).product = A11;
System.reaction(9).propensity = kon1 * A01 * T1;

System.reaction(10).educt = [A10,T1];
System.reaction(10).product = A11;
System.reaction(10).propensity = kon1 * A10 * T1;

System.reaction(11).educt = A11;
System.reaction(11).product = [A11,R];
System.reaction(11).propensity = (r1 + r1) * A11;

System.reaction(12).educt = T1;
System.reaction(12).product = [];
System.reaction(12).propensity = betam1 * T1;

System.reaction(13).educt = [];
System.reaction(13).product = [T1, T1, T1, T1];
System.reaction(13).propensity = beta1;

System.reaction(14).educt = R;
System.reaction(14).product = [];
System.reaction(14).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
