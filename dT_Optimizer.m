% Best Linearization Function
% AUTHOR: Richard Harless
% DATE: 7/14/2026

function [X,Y,MaxDiff,n] = dT_Optimizer(L_Bound,U_Bound,Function,max_err)

%Create Points
n = 10;
too_much_error = 0;
iterations = 0;
increasing = 0;

while too_much_error == 0
    Points = linspace(L_Bound,U_Bound,n);
    relax = 0.01*(U_Bound-L_Bound);
    Max_dT = zeros(n-1,1);
    Graph = plot(Points,Function(Points),'k.','MarkerSize',10);
    hold on
    X1 = linspace(L_Bound,U_Bound,100);
    plot(X1,Function(X1),'b')
    D = 10*ones(n-1,1);
    while max(abs(D)) > 1e-2
        for i = 1:n-1
            Max_dT(i) = dT(Points(i),Points(i+1),Function);
        end
        D = gradient(Max_dT);
        Points(2:n-1) = Points(2:n-1) + transpose(relax*D(2:n-1));
        Graph.XData = Points;
        Graph.YData = Function(Points);
        drawnow
    end
    hold off

    if max(Max_dT) < max_err && n ~= 2
        n = n-1;
        increasing = 0;
    elseif max(Max_dT) > max_err && (iterations == 0 || increasing == 1)
        n = n + 1;
        increasing = 1;
    else
        too_much_error = 1;
    end

    iterations = iterations + 1;
end

X = Points;
Y = Function(Points);
MaxDiff = max(Max_dT);

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
