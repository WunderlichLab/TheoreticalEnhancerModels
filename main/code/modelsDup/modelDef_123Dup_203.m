syms A00 A01 A10 A11 B00 B01 B10 B11 C00 C01 C10 C11 T1 kon1 koff1 r1 myalpha beta1 betam1 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A00; A01; A10; A11; B00; B01; B10; B11; C00; C01; C10; C11; T1; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 10; 150];
System.state.mu0 = [1; 0; 0; 0; 1; 0; 0; 0; 1; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2) + x(3) + x(4)) == 1 && (x(5) + x(6) + x(7) + x(8)) == 1 && (x(9) + x(10) + x(11) + x(12)) == 1);
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

System.reaction(12).educt = B01;
System.reaction(12).product = [B00,T1];
System.reaction(12).propensity = koff1 * B01;

System.reaction(13).educt = B10;
System.reaction(13).product = [B00,T1];
System.reaction(13).propensity = koff1 * B10;

System.reaction(14).educt = [B00,T1];
System.reaction(14).product = B01;
System.reaction(14).propensity = kon1 * B00 * T1;

System.reaction(15).educt = B11;
System.reaction(15).product = [B01,T1];
System.reaction(15).propensity = koff1 * B11;

System.reaction(16).educt = B01;
System.reaction(16).product = [B01,R];
System.reaction(16).propensity = (r1) * B01;

System.reaction(17).educt = [B00,T1];
System.reaction(17).product = B10;
System.reaction(17).propensity = kon1 * B00 * T1;

System.reaction(18).educt = B11;
System.reaction(18).product = [B10,T1];
System.reaction(18).propensity = koff1 * B11;

System.reaction(19).educt = B10;
System.reaction(19).product = [B10,R];
System.reaction(19).propensity = (r1) * B10;

System.reaction(20).educt = [B01,T1];
System.reaction(20).product = B11;
System.reaction(20).propensity = kon1 * B01 * T1;

System.reaction(21).educt = [B10,T1];
System.reaction(21).product = B11;
System.reaction(21).propensity = kon1 * B10 * T1;

System.reaction(22).educt = B11;
System.reaction(22).product = [B11,R];
System.reaction(22).propensity = (r1 + r1) * B11;

System.reaction(23).educt = C01;
System.reaction(23).product = [C00,T1];
System.reaction(23).propensity = koff1 * C01;

System.reaction(24).educt = C10;
System.reaction(24).product = [C00,T1];
System.reaction(24).propensity = koff1 * C10;

System.reaction(25).educt = [C00,T1];
System.reaction(25).product = C01;
System.reaction(25).propensity = kon1 * C00 * T1;

System.reaction(26).educt = C11;
System.reaction(26).product = [C01,T1];
System.reaction(26).propensity = koff1 * C11;

System.reaction(27).educt = C01;
System.reaction(27).product = [C01,R];
System.reaction(27).propensity = (r1) * C01;

System.reaction(28).educt = [C00,T1];
System.reaction(28).product = C10;
System.reaction(28).propensity = kon1 * C00 * T1;

System.reaction(29).educt = C11;
System.reaction(29).product = [C10,T1];
System.reaction(29).propensity = koff1 * C11;

System.reaction(30).educt = C10;
System.reaction(30).product = [C10,R];
System.reaction(30).propensity = (r1) * C10;

System.reaction(31).educt = [C01,T1];
System.reaction(31).product = C11;
System.reaction(31).propensity = kon1 * C01 * T1;

System.reaction(32).educt = [C10,T1];
System.reaction(32).product = C11;
System.reaction(32).propensity = kon1 * C10 * T1;

System.reaction(33).educt = C11;
System.reaction(33).product = [C11,R];
System.reaction(33).propensity = (r1 + r1) * C11;

System.reaction(34).educt = T1;
System.reaction(34).product = [];
System.reaction(34).propensity = betam1 * T1;

System.reaction(35).educt = [];
System.reaction(35).product = [T1, T1, T1, T1];
System.reaction(35).propensity = beta1;

System.reaction(36).educt = R;
System.reaction(36).product = [];
System.reaction(36).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
