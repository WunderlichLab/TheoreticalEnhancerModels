% simulate_enhancer_83.m is the matlab interface to the cvodes mex
%   which simulates the ordinary differential equation and respective
%   sensitivities according to user specifications.
%
% USAGE:
% ======
% [...] = simulate_enhancer_83(tout,theta)
% [...] = simulate_enhancer_83(tout,theta,kappa,options)
% sol = simulate_enhancer_83(...)
% [status,tout,x,y,sx,sy] = simulate_enhancer_83(...)
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
function varargout = simulate_enhancer_83(varargin)

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

options_cvode.ubw = 367;

options_cvode.lbw = 368;

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
rt = [308  309  310  311  312  313  314  315  316  317  318  319  320  321  322  323  324  325  326  327  328  329  330  331  332  333  305   80   85   78   82   83   94   12   23   18   16   13   29   48   59   54   52   49   65  105  116  111  109  106  122  334  335  336   79   95   84   96   91   25   17   30   15   31    1   61   53   66   51   67   37  118  110  123  108  124   73  337  338  339   88   97   89   92    6   32    2   33    3    5   42   68   38   69   39   41   99  125   74  126   75   77  340  341  342   81   98   86   26   19   22   14   24   20   62   55   58   50   60   56  119  112  115  107  117  113  343  344  345   87   93    7   27   21   28   11   34   43   63   57   64   47   70  100  120  114  121  104  127  346  347  348   90   10    9   35    8   36    4   46   45   71   44   72   40  103  102  128  101  129   76  349  350  351  173  178  171  175  176  187  141  152  147  145  142  158  198  209  204  202  199  215  352  353  354  172  188  177  189  184  154  146  159  144  160  130  211  203  216  201  217  166  355  356  357  181  190  182  185  135  161  131  162  132  134  192  218  167  219  168  170  358  359  360  174  191  179  155  148  151  143  153  149  212  205  208  200  210  206  361  362  363  180  186  136  156  150  157  140  163  193  213  207  214  197  220  364  365  366  183  139  138  164  137  165  133  196  195  221  194  222  169  367  368  369  230  235  228  232  233  244  276  287  282  280  277  293  370  371  372  229  245  234  246  241  289  281  294  279  295  223  373  374  375  238  247  239  242  249  296  224  297  225  227  376  377  378  231  248  236  290  283  286  278  288  284  379  380  381  237  243  250  291  285  292  275  298  382  383  384  240  253  252  299  251  300  226  385  386  387  256  261  254  258  259  270  388  389  405  255  271  260  272  267  390  391  392  264  273  265  268  393  394  395  257  274  262  396  397  398  263  269  399  400  401  266  402  403  404  301  302  307  303  306  304];
cw_enhancer_83(sol,tout,theta(1:11),kappa(1:811),options_cvode,plist,pbar,xscale);
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
