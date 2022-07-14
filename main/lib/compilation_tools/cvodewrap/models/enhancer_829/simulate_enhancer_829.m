% simulate_enhancer_829.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%
% USAGE:
% ======
% [...] = simulate_enhancer_829(tout,theta)
% [...] = simulate_enhancer_829(tout,theta,kappa,options)
% sol = simulate_enhancer_829(...)
% [status,tout,x,y,sx,sy] = simulate_enhancer_829(...)
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
function varargout = simulate_enhancer_829(varargin)

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
    kappa(2:419)=0;
   end
else
    kappa=zeros(1,419);
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
if(length(kappa)<419)
    error('provided constant vector is too short');
end

options_cvode.cvodes_atol = 1e-08;
options_cvode.cvodes_rtol = 1e-08;
options_cvode.cvodes_maxsteps = 10000;
options_cvode.sens_ind = 1:11;
options_cvode.nx = 209; % MUST NOT CHANGE THIS VALUE
options_cvode.ny = 2; % MUST NOT CHANGE THIS VALUE
options_cvode.nr = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.ndisc = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.nnz = 4867; % MUST NOT CHANGE THIS VALUE
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

options_cvode.ubw = 206;

options_cvode.lbw = 206;

options_cvode.ss = 0;
options_cvode.qPositiveX = zeros(length(tout),209);

sol.status = 0;
sol.t = tout;
sol.x = zeros(length(tout),209);
sol.y = zeros(length(tout),2);
sol.xdot = zeros(length(tout),209);
sol.J = zeros(209,209);
sol.dydx = zeros(2,209);
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
sol.dxdotdp = zeros(209,options_cvode.np);
plist = options_cvode.sens_ind-1;
if(options_cvode.sensi>0)
    sol.xS = zeros(length(tout),209,length(options_cvode.sens_ind));
    sol.yS = zeros(length(tout),2,length(options_cvode.sens_ind));
    sol.rootS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
    sol.rootvalS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
end
if(max(options_cvode.sens_ind)>11)
    error('Sensitivity index exceeds parameter dimension!')
end
rt = [98   99  100  101  102  103  104  105  106  107  108  109  110  111  112  113  114  115   70   32   42   40   33   43   23   34  116   41   27   35  117   28  118  119    1  120  121  122   24   36  123   25  124  125    2   29  126  127   82  128    3    4  129  130  131  132   37  133   38   47   49   10   39  134  135   83   48   11    7  136  137  138  139   84  140   12    5  141  142   85   86   90    8  143  144   89  145  146  147   26   44   45   18   30  148  149   68   46   19   14  150  151  152  153   20   21  154  155   61   13  156   22  157  158   50  159  160  161   15  162  163   62    6  164   16  165  166   51  167  168  169   52   53  170  171   87  172   54   55  173  174  175  176   31  177  178   77  179   56   57  180  181  182  183   73   78   93   63  184  185   74  186  187   97   79   91    9  188  189   80  190  191  209   88  192   75   81   92   95   96   94   17  193  194   58  195  196  197   64   59  198  199  200  201   60  202  203  204  205   76  206  207  208   65   66   72   67   71   69];
cw_enhancer_829(sol,tout,theta(1:11),kappa(1:419),options_cvode,plist,pbar,xscale);
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
