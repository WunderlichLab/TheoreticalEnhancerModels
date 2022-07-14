% simulate_enhancer_86.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%
% USAGE:
% ======
% [...] = simulate_enhancer_86(tout,theta)
% [...] = simulate_enhancer_86(tout,theta,kappa,options)
% sol = simulate_enhancer_86(...)
% [status,tout,x,y,sx,sy] = simulate_enhancer_86(...)
%
% INPUTS:
% =======
% tout ... 1 dimensional vector of timepoints at which a solution to the ODE is desired
% theta ... 1 dimensional parameter vector of parameters for which sensitivities are desired.
%           this corresponds to the specification in model.sym.p
% kappa ... 1 dimensional parameter vector of parameters for which sensitivities are not desired.
%           this corresponds to the specification in model.sym.k
%           Arbitrary initial conditions can be provided in kappa (see ACME/doc/ACME_doc.pdf for detailed instructions).
% options ... additional options to pass to the cvodes solver. Refer to the cvodes guide for more documentation.
%    .cvodes_atol ... absolute tolerance for the solver. default is specified in the user-provided syms function.
%    .cvodes_rtol ... relative tolerance for the solver. default is specified in the user-provided syms function.
%    .cvodes_maxsteps    ... maximal number of integration steps. default is specified in the user-provided syms function.
%    .tstart    ... start of integration. for all timepoints before this, values will be set to initial value.
%    .sens_ind ... 1 dimensional vector of indexes for which sensitivities must be computed.
%           default value is 1:length(theta).
%    .sx0 ... user-provided sensitivity initialisation. this should be a matrix of dimension [#states x #parameters].
%        default is sensitivity initialisation based on the derivative of the state initialisation.%    .lmm    ... linear multistep method for forward problem.
%        1: Adams-Bashford
%        2: BDF (DEFAULT)
%    .iter    ... iteration method for linear multistep.
%        1: Functional
%        2: Newton (DEFAULT)
%    .linsol   ... linear solver module.
%        direct solvers:
%        1: Dense (DEFAULT)
%        2: Band (not implented)
%        3: LAPACK Dense (not implented)
%        4: LAPACK Band  (not implented)
%        5: Diag (not implented)
%        implicit krylov solvers:
%        6: SPGMR
%        7: SPBCG
%        8: SPTFQMR
%        sparse solvers:
%        9: KLU
%    .stldet   ... flag for stability limit detection. this should be turned on for stiff problems.
%        0: OFF
%        1: ON (DEFAULT)
%    .qPositiveX   ... vector of 0 or 1 of same dimension as state vector. 1 enforces positivity of states.
%    .sensi_meth   ... method for sensitivity computation.
%        1: Forward Sensitivity Analysis (DEFAULT)
%        2: Adjoint Sensitivity Analysis
%    .ism   ... only available for sensi_meth == 1. Method for computation of forward sensitivities.
%        1: Simultaneous (DEFAULT)
%        2: Staggered
%        3: Staggered1
%    .Nd   ... only available for sensi_meth == 2. Number of Interpolation nodes for forward solution. 
%              Default is 1000. 
%    .interpType   ... only available for sensi_meth == 2. Interpolation method for forward solution.
%        1: Hermite (DEFAULT)
%        2: Polynomial
%    .lmmB   ... only available for sensi_meth == 2. linear multistep method for backward problem.
%        1: Adams-Bashford
%        2: BDF (DEFAULT)
%    .iterB   ... only available for sensi_meth == 2. iteration method for linear multistep.
%        1: Functional
%        2: Newton (DEFAULT)
%
% Outputs:
% ========
% sol.status ... flag for status of integration. generally status<0 for failed integration
% sol.tout ... vector at which the solution was computed
% sol.x ... time-resolved state vector
% sol.y ... time-resolved output vector
% sol.sx ... time-resolved state sensitivity vector
% sol.sy ... time-resolved output sensitivity vector
% sol.xdot time-resolved right-hand side of differential equation
% sol.rootval value of root at end of simulation time
% sol.srootval value of root at end of simulation time
% sol.root time of events
% sol.sroot value of root at end of simulation time
function varargout = simulate_enhancer_86(varargin)

% DO NOT CHANGE ANYTHING IN THIS FILE UNLESS YOU ARE VERY SURE ABOUT WHAT YOU ARE DOING
% MANUAL CHANGES TO THIS FILE CAN RESULT IN WRONG SOLUTIONS AND CRASHING OF MATLAB
if(nargin<2)
    error('Not enough input arguments.');
else
    tout=varargin{1};
    phi=varargin{2};
end
if(nargin>=3)
    kappa=varargin{3};
   if(length(kappa)==1)
    kappa(2:505)=0;
   end
else
    kappa=zeros(1,505);
end
if(nargout>1)
    if(nargout>4)
        options_cvode.sensi = 1;
    else
        options_cvode.sensi = 0;
    end
else
    options_cvode.sensi = 1;
end
theta = phi;


if(length(theta)<11)
    error('provided parameter vector is too short');
end
if(length(kappa)<505)
    error('provided constant vector is too short');
end

options_cvode.cvodes_atol = 1e-08;
options_cvode.cvodes_rtol = 1e-08;
options_cvode.cvodes_maxsteps = 10000;
options_cvode.sens_ind = 1:11;
options_cvode.nx = 252; % MUST NOT CHANGE THIS VALUE
options_cvode.ny = 2; % MUST NOT CHANGE THIS VALUE
options_cvode.nr = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.ndisc = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.nnz = 4887; % MUST NOT CHANGE THIS VALUE
options_cvode.tstart = 0;
options_cvode.lmm = 2;
options_cvode.iter = 2;
options_cvode.linsol = 9;
options_cvode.stldet = 1;
options_cvode.Nd = 1000;
options_cvode.interpType = 1;
options_cvode.lmmB = 2;
options_cvode.iterB = 2;
options_cvode.ism = 1;
options_cvode.sensi_meth = 1;

options_cvode.nmaxroot = 100;

options_cvode.ubw = 228;

options_cvode.lbw = 246;

options_cvode.ss = 0;
options_cvode.qPositiveX = zeros(length(tout),252);

sol.status = 0;
sol.t = tout;
sol.x = zeros(length(tout),252);
sol.y = zeros(length(tout),2);
sol.xdot = zeros(length(tout),252);
sol.J = zeros(252,252);
sol.dydx = zeros(2,252);
sol.dydp = zeros(2,11);
sol.root = NaN(options_cvode.nmaxroot,0);
sol.rootval = NaN(options_cvode.nmaxroot,0);
sol.numsteps = zeros(length(tout),1);
sol.numrhsevals = zeros(length(tout),1);
sol.numlinsolvsetups = zeros(length(tout),1);
sol.numerrtestfails = zeros(length(tout),1);
sol.order = zeros(length(tout),1);
sol.numnonlinsolviters = zeros(length(tout),1);
sol.numjacevals = zeros(length(tout),1);
sol.numliniters = zeros(length(tout),1);
sol.numconvfails = zeros(length(tout),1);
sol.numprecevals = zeros(length(tout),1);
sol.numprecsolves = zeros(length(tout),1);

sol.numjtimesevals = zeros(length(tout),1);
sol.numstepsS = zeros(length(tout),1);
sol.numrhsevalsS = zeros(length(tout),1);
sol.numlinsolvsetupsS = zeros(length(tout),1);
sol.numerrtestfailsS = zeros(length(tout),1);
sol.orderS = zeros(length(tout),1);
sol.numnonlinsolvitersS = zeros(length(tout),1);
sol.numjacevalsS = zeros(length(tout),1);
sol.numlinitersS = zeros(length(tout),1);
sol.numconvfailsS = zeros(length(tout),1);
sol.numprecevalsS = zeros(length(tout),1);
sol.numprecsolvesS = zeros(length(tout),1);
sol.numjtimesevalsS = zeros(length(tout),1);

pbar = ones(size(theta));
pbar(pbar==0) = 1;
xscale = [];
if(nargin>=4)
    options_cvode = cw_setdefault(varargin{4},options_cvode);
else
end
if(options_cvode.ss>0)
    options_cvode.sensi = 0;
end
options_cvode.np = length(options_cvode.sens_ind); % MUST NOT CHANGE THIS VALUE
sol.dxdotdp = zeros(252,options_cvode.np);
plist = options_cvode.sens_ind-1;
if(options_cvode.sensi>0)
    sol.xS = zeros(length(tout),252,length(options_cvode.sens_ind));
    sol.yS = zeros(length(tout),2,length(options_cvode.sens_ind));
    sol.rootS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
    sol.rootvalS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
end
if(max(options_cvode.sens_ind)>11)
    error('Sensitivity index exceeds parameter dimension!')
end
rt = [179  180  181  182  183  184  185  186  187  188  189  190  191  192  193  194  195  196  197  198  176   49   51   54   52   47   63   12   25   26    6    7   10   69   82   83   42   43   46  199  200  252   50   53   67   64   55   16   13   15   28   29    8   73   70   72   85   86   44  201  202  203   48   65   66   60   23   18   17   27   30    9   80   75   74   84   87   45  204  205  206   56   57   62   14   24   31   11    1   36   71   81   88   68   37   93  207  208  209   58   61   19   22   32   21    2   33   76   79   89   78   38   90  210  211  212   59   34   20    3   35    5    4   91   77   39   92   41   40  213  214  215  127  129  132  130  125  141  147  160  161  120  121  124  216  217  218  128  131  145  142  133  151  148  150  163  164  122  219  220  221  126  143  144  138  158  153  152  162  165  123  222  223  224  134  135  140  149  159  166  146   94  171  225  226  227  136  139  154  157  167  156   95  168  228  229  230  137  169  155   96  170   98   97  231  232  233  101  103  106  104   99  115  234  235  236  102  105  119  116  107  237  238  239  100  117  118  112  240  241  242  108  109  114  243  244  245  110  113  246  247  248  111  249  250  251  172  173  178  174  177  175];
cw_enhancer_86(sol,tout,theta(1:11),kappa(1:505),options_cvode,plist,pbar,xscale);
sol.x = sol.x(:,rt);
sol.xdot = sol.xdot(:,rt);
if(options_cvode.sensi>0)
    sol.sx = sol.xS(:,rt,:);
    sol.sy = sol.yS;
    sol.sroot = sol.rootS;
    sol.srootval = sol.rootvalS;
elseif(options_cvode.ss>0)
    sol.sx = -sol.J\sol.dxdotdp;
    sol.sy = sol.dydx*sol.sx + sol.dydp;
    sol.sx = sol.sx(rt,:);
end
if(nargout>1)
    varargout{1} = sol.status;
    varargout{2} = sol.t;
    varargout{3} = sol.x;
    varargout{4} = sol.y; % Moments of species
    if(nargout>4)
        varargout{5} = sol.sx;
        varargout{6} = sol.sy;
    end
else
    sol.theta = theta;
    sol.kappa = kappa;
    varargout{1} = sol;
end
end
