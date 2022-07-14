% function plotMM(System,x,y,t,options)
% function plotMM(System,options)
function plotMM(varargin)
if nargin >= 1
    System = varargin{1};
else
    error('At least one input argument is required!')
end
x = System.sol.x;
y = System.sol.y;
t = System.sol.t;

options.save = false;
options.plot_x = true;
options.plot_y = true;
options.plot_xo = false;
options.plot_yo = false;
options.compare = false;
if nargin >= 2
    options = setdefault(varargin{2},options);
end
%% States
if options.plot_x
    % Indices of means and variances
    ind_mean_x = find(sum(System.MM.sym.state.order>=1,2) == 1);
    if System.MM.order >= 2
        ind_var_x = find((System.MM.sym.state.order(:,end-1)== System.MM.sym.state.order(:,end)).*(sum(System.MM.sym.state.order~=0,2)==2));
    end
    
    mean_x = x(:,ind_mean_x);
    if System.MM.order >= 2
        % var_x = x(:,ind_var_x);
        % Alvaro edit
        var_x = sqrt(x(:,ind_var_x));
    else
        var_x = [];
    end
    
    options_x = options;
    options_x.ylabel = System.state.name;
    if options.compare
        options_x.fig_title = {options.fig_title_x{1}};
        options_x.fh = options.fhx{1};
        plotLine(mean_x,t,options_x);
        options_x.fig_title = {options.fig_title_x{2}};
        options_x.fh = options.fhx{2};
        plotLine(var_x,t,options_x);
    else
        options_x.fig_title = {'MM states'};
        plotInterval(mean_x,var_x,t,options_x);
    end
end
%% Outputs
if isfield(System,'output')
    if options.plot_y
        % Indices of means and variances
        ind_mean_y = find(sum(System.MM.sym.output.order>=1,2) == 1);
        if System.MM.output_order >= 2
            ind_var_y = find((sum(System.MM.sym.output.order~=0,2)==2).*(System.MM.sym.output.order(:,end-1)== System.MM.sym.output.order(:,end)));
        end
        
        mean_y = y(:,ind_mean_y);
        if System.MM.output_order >= 2
            var_y = sqrt(y(:,ind_var_y));
        else
            var_y = [];
        end
        
        options_y = options;
        options_y.ylabel = System.output.name;
        if options.compare
            options_y.fig_title = {options.fig_title_y{1}};
            options_y.fh = options.fhy{1};
            plotLine(mean_y,t,options_y);
            options_y.fig_title = {options.fig_title_y{2}};
            options_y.fh = options.fhy{2};
            plotLine(var_y,t,options_y);
        else
            options_y.fig_title = {'MM outputs'};
            plotInterval(mean_y,var_y,t,options_y);
        end
    end
end
%% Higher-order moments
if System.MM.order >= 2
    if options.plot_xo
        if isfield(options,'state_order')
            xo = options.state_order;
            options_xo = options;
            ls_species = [];
            for i = 1:System.state.number-1
                ls_species = [ls_species,num2str(i),': ',System.state.name{i},',   '];
            end
            ls_species = [ls_species,num2str(System.state.number),': ',System.state.name{end}];
            options_xo.fs = 12;
            for i = 2:xo
                options_xo.fig_title = {['Species moments of order ',num2str(i)],ls_species};
                ind_xo = find((sum(System.MM.sym.state.order~=0,2)==i));
                options_xo.ylabel = System.MM.sym.state.moments(ind_xo);
                plotLine(x(:,ind_xo),t,options_xo)
            end
        end
    end
end
if isfield(System,'output')
    if System.MM.output_order >= 2
        if options.plot_yo
            if isfield(options,'output_order')
                yo = options.output_order;
                options_yo = options;
                ls_outputs = [];
                for i = 1:System.output.number-1
                    ls_outputs = [ls_outputs,num2str(i),': ',System.output.name{i},',   '];
                end
                ls_outputs = [ls_outputs,num2str(System.output.number),': ',System.output.name{end}];
                options_yo.fs = 12;
                for i = 2:yo
                    options_yo.fig_title = {['Outputs moments of order ',num2str(i)],ls_outputs};
                    ind_yo = find((sum(System.MM.sym.output.order~=0,2)==i));
                    options_yo.ylabel = System.MM.sym.output.moments(ind_yo);
                    plotLine(y(:,ind_yo),t,options_yo)
                end
            end
        end
    end
end