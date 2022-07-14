syms A0 A1 B000 B002 B020 B022 B200 B202 B220 B222 T1 T2 kon1 koff1 r1 kon2 koff2 r2 myalpha beta1 betam1 beta2 betam2 Omega time R

System.time = time;
System.compartments = {'cell'};
System.volumes = [Omega];
System.state.variable = [A0; A1; B000; B002; B020; B022; B200; B202; B220; B222; T1; T2; R];
System.state.compartment = {'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'; 'cell'};
System.state.type = {'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'stochastic';'moment'};
System.state.xmin = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.xmax = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 10; 10; 150];
System.state.mu0 = [1; 0; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
System.state.C0 = zeros(length(System.state.variable)*(length(System.state.variable)+1)/2,1);
System.state.constraint = @(x) ((x(1) + x(2)) == 1 && (x(3) + x(4) + x(5) + x(6) + x(7) + x(8) + x(9) + x(10)) == 1);
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

System.reaction(4).educt = B002;
System.reaction(4).product = [B000,T2];
System.reaction(4).propensity = koff2 * B002;

System.reaction(5).educt = B020;
System.reaction(5).product = [B000,T2];
System.reaction(5).propensity = koff2 * B020;

System.reaction(6).educt = B200;
System.reaction(6).product = [B000,T2];
System.reaction(6).propensity = koff2 * B200;

System.reaction(7).educt = [B000,T2];
System.reaction(7).product = B002;
System.reaction(7).propensity = kon2 * B000 * T2;

System.reaction(8).educt = B022;
System.reaction(8).product = [B002,T2];
System.reaction(8).propensity = koff2 * B022;

System.reaction(9).educt = B202;
System.reaction(9).product = [B002,T2];
System.reaction(9).propensity = koff2 * B202;

System.reaction(10).educt = B002;
System.reaction(10).product = [B002,R];
System.reaction(10).propensity = (r2) * B002;

System.reaction(11).educt = [B000,T2];
System.reaction(11).product = B020;
System.reaction(11).propensity = kon2 * B000 * T2;

System.reaction(12).educt = B022;
System.reaction(12).product = [B020,T2];
System.reaction(12).propensity = koff2 * B022;

System.reaction(13).educt = B220;
System.reaction(13).product = [B020,T2];
System.reaction(13).propensity = koff2 * B220;

System.reaction(14).educt = B020;
System.reaction(14).product = [B020,R];
System.reaction(14).propensity = (r2) * B020;

System.reaction(15).educt = [B000,T2];
System.reaction(15).product = B200;
System.reaction(15).propensity = kon2 * B000 * T2;

System.reaction(16).educt = B202;
System.reaction(16).product = [B200,T2];
System.reaction(16).propensity = koff2 * B202;

System.reaction(17).educt = B220;
System.reaction(17).product = [B200,T2];
System.reaction(17).propensity = koff2 * B220;

System.reaction(18).educt = B200;
System.reaction(18).product = [B200,R];
System.reaction(18).propensity = (r2) * B200;

System.reaction(19).educt = [B002,T2];
System.reaction(19).product = B022;
System.reaction(19).propensity = kon2 * B002 * T2;

System.reaction(20).educt = [B020,T2];
System.reaction(20).product = B022;
System.reaction(20).propensity = kon2 * B020 * T2;

System.reaction(21).educt = B222;
System.reaction(21).product = [B022,T2];
System.reaction(21).propensity = koff2 * B222;

System.reaction(22).educt = B022;
System.reaction(22).product = [B022,R];
System.reaction(22).propensity = (r2 + r2) * B022;

System.reaction(23).educt = [B002,T2];
System.reaction(23).product = B202;
System.reaction(23).propensity = kon2 * B002 * T2;

System.reaction(24).educt = [B200,T2];
System.reaction(24).product = B202;
System.reaction(24).propensity = kon2 * B200 * T2;

System.reaction(25).educt = B222;
System.reaction(25).product = [B202,T2];
System.reaction(25).propensity = koff2 * B222;

System.reaction(26).educt = B202;
System.reaction(26).product = [B202,R];
System.reaction(26).propensity = (r2 + r2) * B202;

System.reaction(27).educt = [B020,T2];
System.reaction(27).product = B220;
System.reaction(27).propensity = kon2 * B020 * T2;

System.reaction(28).educt = [B200,T2];
System.reaction(28).product = B220;
System.reaction(28).propensity = kon2 * B200 * T2;

System.reaction(29).educt = B222;
System.reaction(29).product = [B220,T2];
System.reaction(29).propensity = koff2 * B222;

System.reaction(30).educt = B220;
System.reaction(30).product = [B220,R];
System.reaction(30).propensity = (r2 + r2) * B220;

System.reaction(31).educt = [B022,T2];
System.reaction(31).product = B222;
System.reaction(31).propensity = kon2 * B022 * T2;

System.reaction(32).educt = [B202,T2];
System.reaction(32).product = B222;
System.reaction(32).propensity = kon2 * B202 * T2;

System.reaction(33).educt = [B220,T2];
System.reaction(33).product = B222;
System.reaction(33).propensity = kon2 * B220 * T2;

System.reaction(34).educt = B222;
System.reaction(34).product = [B222,R];
System.reaction(34).propensity = (r2 + r2 + r2) * B222;

System.reaction(35).educt = T1;
System.reaction(35).product = [];
System.reaction(35).propensity = betam1 * T1;

System.reaction(36).educt = [];
System.reaction(36).product = [T1, T1, T1, T1];
System.reaction(36).propensity = beta1;

System.reaction(37).educt = T2;
System.reaction(37).product = [];
System.reaction(37).propensity = betam2 * T2;

System.reaction(38).educt = [];
System.reaction(38).product = [T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2, T2];
System.reaction(38).propensity = beta2;

System.reaction(39).educt = R;
System.reaction(39).product = [];
System.reaction(39).propensity = myalpha * R;

System.output.variable = [R];
System.output.function = [R];
