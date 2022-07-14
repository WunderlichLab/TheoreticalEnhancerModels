syms A0 A1 B00 B02 B10 B12 T1 T2 kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0; A1; B00; B02; B10; B12; T1; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 10; 10; 150];
System.state.mu0 = [1; 0; 1; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2)) == 1 && (x(3) + x(4) + x(5) + x(6)) == 1);
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

System.reaction(4).educt = B02;
System.reaction(4).product = [B00,T2];
System.reaction(4).propensity = koff2 * B02;

System.reaction(5).educt = B10;
System.reaction(5).product = [B00,T1];
System.reaction(5).propensity = koff1 * B10;

System.reaction(6).educt = [B00,T2];
System.reaction(6).product = B02;
System.reaction(6).propensity = kon2 * B00 * T2;

System.reaction(7).educt = B12;
System.reaction(7).product = [B02,T1];
System.reaction(7).propensity = koff1 * B12;

System.reaction(8).educt = B02;
System.reaction(8).product = [B02,R];
System.reaction(8).propensity = (r2) * B02;

System.reaction(9).educt = [B00,T1];
System.reaction(9).product = B10;
System.reaction(9).propensity = kon1 * B00 * T1;

System.reaction(10).educt = B12;
System.reaction(10).product = [B10,T2];
System.reaction(10).propensity = koff2 * B12;

System.reaction(11).educt = B10;
System.reaction(11).product = [B10,R];
System.reaction(11).propensity = (r1) * B10;

System.reaction(12).educt = [B02,T1];
System.reaction(12).product = B12;
System.reaction(12).propensity = kon1 * B02 * T1;

System.reaction(13).educt = [B10,T2];
System.reaction(13).product = B12;
System.reaction(13).propensity = kon2 * B10 * T2;

System.reaction(14).educt = B12;
System.reaction(14).product = [B12,R];
System.reaction(14).propensity = (r1 + r2) * B12;

System.reaction(15).educt = T1;
System.reaction(15).product = [];
System.reaction(15).propensity = betam1 * T1;

System.reaction(16).educt = [];
System.reaction(16).product = [T1, T1, T1, T1];
System.reaction(16).propensity = beta1;

System.reaction(17).educt = T2;
System.reaction(17).product = [];
System.reaction(17).propensity = betam2 * T2;

System.reaction(18).educt = [];
System.reaction(18).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(18).propensity = beta2;

System.reaction(19).educt = R;
System.reaction(19).product = [];
System.reaction(19).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
