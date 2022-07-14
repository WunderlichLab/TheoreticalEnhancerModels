% simulate_enhancer_87.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%
% USAGE:
% ======
% [...] = simulate_enhancer_87(tout,theta)
% [...] = simulate_enhancer_87(tout,theta,kappa,options)
% sol = simulate_enhancer_87(...)
% [status,tout,x,y,sx,sy] = simulate_enhancer_87(...)
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
function varargout = simulate_enhancer_87(varargin)

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
    kappa(2:811)=0;
   end
else
    kappa=zeros(1,811);
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
if(length(kappa)<811)
    error('provided constant vector is too short');
end

options_cvode.cvodes_atol = 1e-08;
options_cvode.cvodes_rtol = 1e-08;
options_cvode.cvodes_maxsteps = 10000;
options_cvode.sens_ind = 1:11;
options_cvode.nx = 405; % MUST NOT CHANGE THIS VALUE
options_cvode.ny = 2; % MUST NOT CHANGE THIS VALUE
options_cvode.nr = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.ndisc = 0; % MUST NOT CHANGE THIS VALUE
options_cvode.nnz = 8238; % MUST NOT CHANGE THIS VALUE
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

options_cvode.ubw = 363;

options_cvode.lbw = 399;

options_cvode.ss = 0;
options_cvode.qPositiveX = zeros(length(tout),405);

sol.status = 0;
sol.t = tout;
sol.x = zeros(length(tout),405);
sol.y = zeros(length(tout),2);
sol.xdot = zeros(length(tout),405);
sol.J = zeros(405,405);
sol.dydx = zeros(2,405);
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
sol.dxdotdp = zeros(405,options_cvode.np);
plist = options_cvode.sens_ind-1;
if(options_cvode.sensi>0)
    sol.xS = zeros(length(tout),405,length(options_cvode.sens_ind));
    sol.yS = zeros(length(tout),2,length(options_cvode.sens_ind));
    sol.rootS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
    sol.rootvalS =  NaN(options_cvode.nmaxroot,0,length(options_cvode.sens_ind));
end
if(max(options_cvode.sens_ind)>11)
    error('Sensitivity index exceeds parameter dimension!')
end
rt = [308  309  310  311  312  313  314  315  316  317  318  319  320  321  322  323  324  325  326  327  328  329  330  331  332  333  305   85   87   90   88   83   99   12   25   26    6    7   10   48   61   62   42   43   46  105  118  119   78   79   82  334  335  405   86   89  103  100   91   16   13   15   28   29    8   52   49   51   64   65   44  109  106  108  121  122   80  336  337  338   84  101  102   96   23   18   17   27   30    9   59   54   53   63   66   45  116  111  110  120  123   81  339  340  341   92   93   98   14   24   31   11    1   36   50   60   67   47   37   72  107  117  124  104   73  129  342  343  344   94   97   19   22   32   21    2   33   55   58   68   57   38   69  112  115  125  114   74  126  345  346  347   95   34   20    3   35    5    4   70   56   39   71   41   40  127  113   75  128   77   76  348  349  350  178  180  183  181  176  192  141  154  155  135  136  139  198  211  212  171  172  175  351  352  353  179  182  196  193  184  145  142  144  157  158  137  202  199  201  214  215  173  354  355  356  177  194  195  189  152  147  146  156  159  138  209  204  203  213  216  174  357  358  359  185  186  191  143  153  160  140  130  165  200  210  217  197  166  222  360  361  362  187  190  148  151  161  150  131  162  205  208  218  207  167  219  363  364  365  188  163  149  132  164  134  133  220  206  168  221  170  169  366  367  368  256  258  261  259  254  270  276  289  290  249  250  253  369  370  371  257  260  274  271  262  280  277  279  292  293  251  372  373  374  255  272  273  267  287  282  281  291  294  252  375  376  377  263  264  269  278  288  295  275  223  300  378  379  380  265  268  283  286  296  285  224  297  381  382  383  266  298  284  225  299  227  226  384  385  386  230  232  235  233  228  244  387  388  389  231  234  248  245  236  390  391  392  229  246  247  241  393  394  395  237  238  243  396  397  398  239  242  399  400  401  240  402  403  404  301  302  307  303  306  304];
cw_enhancer_87(sol,tout,theta(1:11),kappa(1:811),options_cvode,plist,pbar,xscale);
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
