% Best Linearization Function
% AUTHOR: Richard Harless
% DATE: 7/14/2026

function [X,Y,MaxDiff,n_points] = dT_Optimizer(L_Bound,U_Bound,max_err,Function)

arguments
    L_Bound (1,1) double % Lower Bound of Input Values
    U_Bound (1,1) double % Upper Bound of Input Values
    max_err (1,1) double % Maximum allowed percent error
end
arguments (Repeating)
    Function (1,1) function_handle % Function
end

%Create Points
n_points = 10;
n_functions = numel(Function);
too_much_error = 0;
iterations = 0;
increasing = 0;

while too_much_error == 0
    Points = linspace(L_Bound,U_Bound,n_points);
    relax = 0.03*(U_Bound-L_Bound);

    Max_dT = 10*ones(n_points-1,n_functions);

    % Graph = plot(Points,Function(Points),'ko','MarkerSize',6,'MarkerFaceColor','k');
    % hold on
    % X1 = linspace(L_Bound,U_Bound,100);
    % plot(X1,Function(X1),'b')

    D = 10*ones(n_points-1,n_functions);
    Delta = 10*ones(n_points-1,1);
    while max(abs(Delta),[],'all') > 1e-3 && max(Max_dT,[],'all') > max_err
        for i = 1:n_points-1
            for j = 1:n_functions
            Max_dT(i,j) = dT(Points(i),Points(i+1),Function{j});
            end
        end
        for j = 1:n_functions
            D(:,j) = gradient(Max_dT(:,j));
        end
        Points(2:n_points-1) = Points(2:n_points-1) + transpose(relax*sum(D(2:n_points-1,:),2));
        Delta = transpose(relax*sum(D(2:n_points-1,:),2));

        % Graph.XData = Points;
        % Graph.YData = Function(Points);
        % drawnow
    end
    % hold off

    if max(Max_dT,[],'all') < max_err && n_points ~= 2
        n_points = n_points-1;
        increasing = 0;
    elseif max(Max_dT,[],'all') > max_err && (iterations == 0 || increasing == 1)
        n_points = n_points + 1;
        increasing = 1;
    else
        too_much_error = 1;
    end

    iterations = iterations + 1;
end

X = Points;
Y = zeros(n_points,n_functions);
for i = 1:n_functions
    Y(:,i) = Function{i}(Points);
end
MaxDiff = max(Max_dT,[],'all');

    function max_dT_percent = dT(A,B,Function)

        y_A = Function(A);
        y_B = Function(B);

        m = (y_B - y_A)/(B - A);
        linear = @(x1) m*(x1 - A) + y_A;
        percent_diff_inv = @(x2) -1*(abs(Function(x2) - linear(x2)))/Function(x2);
        [~,max_dT_percent_inv] = fminbnd(percent_diff_inv,A,B);
        max_dT_percent = -100*max_dT_percent_inv;
    end

end
