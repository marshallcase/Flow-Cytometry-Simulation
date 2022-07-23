function frac_bound = getFracBound(kon1,koff1,kon2,koff2,P1,P2)
    %initial values
    L_0 = 1e10/(6.022e23*25e-6); %M
    LP1_0 = 0; %M
    LP2_0 = 0; %M

    %start time
    t_0 = 0; %s
    t_end = 3600*25; %s

    %number of time steps
    n=5001;
    
    %time
    time = linspace(t_0,t_end,n);
    
    %get numerical solutions
    [x,c] = ode15s(@(x,c)coupled(x,c,kon1,koff1,kon2,koff2,P1,P2),...
    time,[L_0,LP1_0,LP2_0]);

    frac_bound = c(end,2)/L_0;
end

function dcdt = coupled(x,c,kon1,koff1,kon2,koff2,P1,P2)
%differential equations governing peptide behavior
dcdt = [-kon1*c(1)*P1+koff1*c(2)-kon2*c(1)*P2+koff2*c(3);...
            kon1*c(1)*P1-koff1*c(2);...
            kon2*c(1)*P2-koff2*c(3)];
end