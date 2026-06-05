function cardiac_curve = Make_Cardiac_Curve(period, tot_samples)
    %{
    Create a cardiac curve that matches literature and downsample for our framerate
    2/18/26 Kyler Howard

    param: period         - Number of samples each cardiac curve should take
    param: tot_samples    - Number of samples to fill with cardiac curves

    return: cardiac_curve - Repeated cardiac curve normalized from 0-1
    %}

    % Endpoints in time scale
    a = 0;
    b = 0.06;
    c = 0.25;
    d = 0.31;
    e = 0.8;
    
    % Create original signal on cardiac scale x/y
    Np = 1000;
    x1 = linspace(a,e,Np);
    
    % Exponential coefficients fitting to (0.06,120), (0.2, 53), (0.25, 50)
    a2 = 51.1928;
    b2 = 347.8166;
    c2 = -26.29949;
    
    % Logrithmic coefficients fitting to (0.31, 50), (0.51, 102), (0.80, 120)
    a4 = 168.10094202;
    b4 = -36.09097339;
    
    % Create exponential and logrithmic curves
    y2 = a2 + b2*exp(c2*x1(x1>=b & x1<=c));
    y4 = a4 + b4 ./ x1(x1>=d & x1<=e);
    
    % Fit the flat points to exactly connect between fitted curves
    m1 = (y2(1)-y4(end)) / (b-a);
    b1 = y2(1) - m1*b;
    m3 = (y4(1)-y2(end)) / (d-c);
    b3 = y4(1) - m3*d;
    
    y1 = m1*x1(x1>=a & x1<=b) + b1;
    y3 = m3*x1(x1>=c & x1<=d) + b3;
    
    % Create cardiac curve on cardiac scales
    yold = zeros(size(x1));
    yold(x1>=a & x1<=b) = y1;
    yold(x1>=b & x1<=c) = y2;
    yold(x1>=c & x1<=d) = y3;
    yold(x1>=d & x1<=e) = y4;
    clear y1 y2 y3 y4
    
    % Original domain length
    L = e - a;
    
    % Normalized domain
    u = linspace(0,1,Np);
    x = a + L*u;
    
    % Normalized breakpoints
    ub = b/L;
    uc = c/L;
    ud = d/L;
    
    % Preallocate
    y = zeros(size(u));
    
    % Create segments in u domain
    idx1 = (u >= 0 & u <= ub);
    idx2 = (u > ub & u <= uc);
    idx3 = (u > uc & u <= ud);
    idx4 = (u > ud & u <= 1);

    y(idx1) = m1*x(idx1) + b1;
    y(idx2) = a2 + b2*exp(c2*x(idx2));
    y(idx3) = m3*x(idx3) + b3;
    y(idx4) = a4 + b4 ./ x(idx4);
    
    % Normalize from 0-1
    y = (y - min(y)) / (max(y) - min(y));
    
    % Random shift of the curve & domain
    k         = randi(Np);
    u_shift   = (k-1)/Np;
    u_shifted = mod(u + u_shift, 1);
    [~, iq]   = sort(u_shifted);
    y_shifted = y(iq);
    
    % figure(1);
    %     clf
    %     subplot(2,1,1)
    %         plot(x1, yold)
    %     subplot(2,1,2)
    %         hold on
    %         plot(u, y)
    %         scatter(u_shifted, y_shifted,"filled")
    %         legend("Regular", "Shifted")
    %         ylim([-0.1, 1.1])
    
    % Downsample with interopolation
    u_sub = linspace(0,1,period);
    y_sub = interp1(u_shifted(1:Np-1), y_shifted(1:Np-1), u_sub, "linear", "extrap");
    
    % Repeat the signal to fill out breaths
    n_periods = floor(tot_samples / (period-1)); 
    
    % Repeat full periods 
    cardiac_curve = repmat(y_sub(2:end), 1, n_periods); 
    
    % Fill remaining samples 
    remaining = tot_samples - length(cardiac_curve); 
    if remaining > 0 
        cardiac_curve = [cardiac_curve, y_sub(2:remaining+1)]; 
    end
    
    % figure(2)
    %     clf
    %     hold on
    %     plot(1:frames, y_rep)
    %     plot([14,14], [0,1], ":r")
    %     ylim([-0.1, 1.1])
end
