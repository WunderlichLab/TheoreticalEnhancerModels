% function [X,Y,XRed,SX,SY] = simMM2_reduced(t,theta,x0) 
function varargout = simMM2_reduced(varargin) 

t = varargin{1};
theta = varargin{2};
x0 = [];
if nargin>2
    x0 = varargin{3};
end
data.theta = theta;
% Set solver options 
options_CVode = CVodeSetOptions('RelTol',1e-6,...
                                'AbsTol',1e-6,...
                                'MaxNumSteps',10^6,...
                                'JacobianFn',@jacfn,...
                                'Userdata',data);
options_CVodes = CVodeSensSetOptions('method','Simultaneous',...
                                     'ErrControl',true,...
                                     'ParamScales',1:length(theta));

% Initial conditions
if isempty(x0)
    x0 = x0fun(theta);
end
if nargout >= 4
    sx0 = sx0fun(theta);
end

% Initialization of CVode
CVodeInit(@rhs,'BDF','Newton',0,x0,options_CVode);
if nargout >= 4
    CVodeSensInit(length(theta),@rhsS,sx0,options_CVodes);
end

% Simulation
if nargout <= 3
    [status,~,x] = CVode(t(2:end),'Normal');
    X = [x0';x'];
    Y = rhsO(t,X,theta);
elseif nargout >=4
    [status,~,x,sx] = CVode(t(2:end),'Normal');
    X = [x0';x'];
    Y = rhsO(t,X,theta);
    SX = zeros(length(t),length(x0),length(theta));
    SX(1,:,:) = sx0;
    SX(2:end,:,:) = permute(sx,[3,1,2]);
    SY = rhsOS(t,X,SX,Y,theta);
end

% Evaluate reduced covariances
if nargout >= 3
XRed = EvalRedCov(X);
end
% Free memory
CVodeFree;

% Assign output
varargout{1} = X;
if nargout >= 2
    varargout{2} = Y;
end
if nargout >= 3
    varargout{3} = XRed;
end
if nargout >= 4
    varargout{4} = SX;
end
if nargout >= 5
    varargout{5} = SY;
end
if nargout >= 6
    error('Too many output arguments.');
end


%% RIGHT-HAND SIDE
function [dxdt,flag,new_data] = rhs(t,x,data) 

theta = data.theta;
dxdt = [-theta(1)*x(1);...
        theta(1)*x(1)-theta(2)*x(2);...
        theta(2)*x(2)+theta(4)*x(4);...
        -theta(4)*x(4);...
        theta(1)*x(1)-2*x(5)*theta(1);...
        x(5)*theta(1)-x(6)*theta(1)-x(6)*theta(2)-theta(1)*x(1);...
        2*x(6)*theta(1)-2*x(7)*theta(2)+theta(1)*x(1)+theta(2)*x(2);...
        x(7)*theta(2)-x(8)*theta(2)-theta(2)*x(2)+(x(6)*x(8)*theta(1))/(x(7)+1/10000000000)+(x(6)^2*x(8)*x(10)*theta(4))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000));...
        2*x(8)*theta(2)+2*x(10)*theta(4)+theta(2)*x(2)+theta(4)*x(4);...
        x(11)*theta(4)-x(10)*theta(4)-theta(4)*x(4)+(x(6)^2*x(8)*x(10)*theta(2))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000));...
        theta(4)*x(4)-2*x(11)*theta(4)];

flag = 0;
new_data = [];

%% RIGHT-HAND SIDE OF SENSITIVITIES
function [dsxdt,flag,new_data] = rhsS(t,x,dxdt,sx,data) 

theta = data.theta;
J = jacfn(t,x,dxdt,data);
dfdtheta = [-x(1),0,0,0;...
              x(1),-x(2),0,0;...
              0,x(2),0,x(4);...
              0,0,0,-x(4);...
              x(1)-2*x(5),0,0,0;...
              x(5)-x(6)-x(1),-x(6),0,0;...
              2*x(6)+x(1),x(2)-2*x(7),0,0;...
              (x(6)*x(8))/(x(7)+1/10000000000),x(7)-x(8)-x(2),0,(x(6)^2*x(8)*x(10))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000));...
              0,2*x(8)+x(2),0,2*x(10)+x(4);...
              0,(x(6)^2*x(8)*x(10))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)),0,x(11)-x(10)-x(4);...
              0,0,0,x(4)-2*x(11)];

dsxdt = J*sx + dfdtheta;

flag = 0;
new_data = [];

%% JACOBIAN
function [J,flag,new_data] = jacfn(t,x,dxdt,data) 

theta = data.theta;
J = [-theta(1),0,0,0,0,0,0,0,0,0,0;...
     theta(1),-theta(2),0,0,0,0,0,0,0,0,0;...
     0,theta(2),0,theta(4),0,0,0,0,0,0,0;...
     0,0,0,-theta(4),0,0,0,0,0,0,0;...
     theta(1),0,0,0,-2*theta(1),0,0,0,0,0,0;...
     -theta(1),0,0,0,theta(1),-theta(1)-theta(2),0,0,0,0,0;...
     theta(1),theta(2),0,0,0,2*theta(1),-2*theta(2),0,0,0,0;...
     0,-theta(2),0,0,-(x(6)^2*x(8)*x(10)*theta(4))/((x(5)+1/10000000000)^2*(x(7)+1/10000000000)*(x(9)+1/10000000000)),(x(8)*theta(1))/(x(7)+1/10000000000)+(2*x(6)*x(8)*x(10)*theta(4))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)),theta(2)-(x(6)*x(8)*theta(1))/(x(7)+1/10000000000)^2-(x(6)^2*x(8)*x(10)*theta(4))/((x(5)+1/10000000000)*(x(7)+1/10000000000)^2*(x(9)+1/10000000000)),(x(6)*theta(1))/(x(7)+1/10000000000)-theta(2)+(x(6)^2*x(10)*theta(4))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)),-(x(6)^2*x(8)*x(10)*theta(4))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)^2),(x(6)^2*x(8)*theta(4))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)),0;...
     0,theta(2),0,theta(4),0,0,0,2*theta(2),0,2*theta(4),0;...
     0,0,0,-theta(4),-(x(6)^2*x(8)*x(10)*theta(2))/((x(5)+1/10000000000)^2*(x(7)+1/10000000000)*(x(9)+1/10000000000)),(2*x(6)*x(8)*x(10)*theta(2))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)),-(x(6)^2*x(8)*x(10)*theta(2))/((x(5)+1/10000000000)*(x(7)+1/10000000000)^2*(x(9)+1/10000000000)),(x(6)^2*x(10)*theta(2))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)),-(x(6)^2*x(8)*x(10)*theta(2))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000)^2),(x(6)^2*x(8)*theta(2))/((x(5)+1/10000000000)*(x(7)+1/10000000000)*(x(9)+1/10000000000))-theta(4),theta(4);...
     0,0,0,theta(4),0,0,0,0,0,0,-2*theta(4)];
flag = 0;
new_data = [];

%% OUTPUT MAP
function y = rhsO(t,x,theta) 

y = [];


%% OUTPUT MAP OF SENSITIVITIES
function sy = rhsOS(t,x,sx,y,theta) 

sy = zeros(length(t),size(y,2),length(theta));
for k = 1:length(t)
    dHdx = [];
    dHdtheta = [];
    sy(k,:,:) = dHdx*squeeze(sx(k,:,:)) + dHdtheta;
end


%% INITIAL CONDITIONS FOR STATE
function x0 = x0fun(theta) 

x0 = [20;...
      0;...
      0;...
      0;...
      0;...
      0;...
      0;...
      0;...
      0;...
      0;...
      0];


%% INITIAL CONDITIONS FOR STATE SENSITIVITY
function sx0 = sx0fun(theta) 

sx0 = [0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0;...
       0,0,0,0];


%% EVALUAED REDUCED COVARIANCES
function xred = EvalRedCov(x) 

xred =[(x(:,6).*x(:,8))./(x(:,7)+1./10000000000),(x(:,6).*x(:,8).*x(:,10))./((x(:,7)+1./10000000000).*(x(:,9)+1./10000000000)),(x(:,6).^2.*x(:,8).*x(:,10))./((x(:,5)+1./10000000000).*(x(:,7)+1./10000000000).*(x(:,9)+1./10000000000))];


