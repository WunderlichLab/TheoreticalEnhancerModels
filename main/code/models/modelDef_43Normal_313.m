syms A0 A1 B0 B1 C00 C02 C10 C12 T1 T2 kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0; A1; B0; B1; C00; C02; C10; C12; T1; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 10; 10; 150];
System.state.mu0 = [1; 0; 1; 0; 1; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2)) == 1 && (x(3) + x(4)) == 1 && (x(5) + x(6) + x(7) + x(8)) == 1);
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

System.reaction(7).educt = C02;
System.reaction(7).product = [C00,T2];
System.reaction(7).propensity = koff2 * C02;

System.reaction(8).educt = C10;
System.reaction(8).product = [C00,T1];
System.reaction(8).propensity = koff1 * C10;

System.reaction(9).educt = [C00,T2];
System.reaction(9).product = C02;
System.reaction(9).propensity = kon2 * C00 * T2;

System.reaction(10).educt = C12;
System.reaction(10).product = [C02,T1];
System.reaction(10).propensity = koff1 * C12;

System.reaction(11).educt = C02;
System.reaction(11).product = [C02,R];
System.reaction(11).propensity = (r2) * C02;

System.reaction(12).educt = [C00,T1];
System.reaction(12).product = C10;
System.reaction(12).propensity = kon1 * C00 * T1;

System.reaction(13).educt = C12;
System.reaction(13).product = [C10,T2];
System.reaction(13).propensity = koff2 * C12;

System.reaction(14).educt = C10;
System.reaction(14).product = [C10,R];
System.reaction(14).propensity = (r1) * C10;

System.reaction(15).educt = [C02,T1];
System.reaction(15).product = C12;
System.reaction(15).propensity = kon1 * C02 * T1;

System.reaction(16).educt = [C10,T2];
System.reaction(16).product = C12;
System.reaction(16).propensity = kon2 * C10 * T2;

System.reaction(17).educt = C12;
System.reaction(17).product = [C12,R];
System.reaction(17).propensity = (r1 + r2) * C12;

System.reaction(18).educt = T1;
System.reaction(18).product = [];
System.reaction(18).propensity = betam1 * T1;

System.reaction(19).educt = [];
System.reaction(19).product = [T1, T1, T1, T1];
System.reaction(19).propensity = beta1;

System.reaction(20).educt = T2;
System.reaction(20).product = [];
System.reaction(20).propensity = betam2 * T2;

System.reaction(21).educt = [];
System.reaction(21).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(21).propensity = beta2;

System.reaction(22).educt = R;
System.reaction(22).product = [];
System.reaction(22).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
