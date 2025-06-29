function gauss_seidel_14bus()
    % IEEE 14-Bus System Power Flow Analysis using Gauss-Seidel Method
    
    bus_data = [
        1   1.060   0   114.17   -16.9   0       0       0       10;
        2   1.045   0   40.00    0      21.7     12.7    -42.0   50.0;
        3   1.010   0   0        0      94.2     19.1    23.4    40.0;
        4   1.000   0   0        0      47.8     -3.9    0       0;
        5   1.000   0   0        0      7.6      1.6     0       0;
        6   1.000   0   0        0      11.2     7.5     0       0;
        7   1.000   0   0        0      0        0       0       0;
        8   1.000   0   0        0      0        0       0       0;
        9   1.000   0   0        0      29.5     16.6    0       0;
        10  1.000   0   0        0      9.0      5.8     0       0;
        11  1.000   0   0        0      3.5      1.8     0       0;
        12  1.000   0   0        0      6.1      1.6     0       0;
        13  1.000   0   0        0      13.8     5.8     0       0;
        14  1.000   0   0        0      14.9     5.0     0       0;
    ];
    
    % Line data from Table A.1
    line_data = [
        1   1   2   0.01938   0.05917   0.02640;
        2   1   5   0.05403   0.22304   0.02190;
        3   2   3   0.04699   0.19797   0.01870;
        4   2   4   0.05811   0.17632   0.02460;
        5   2   5   0.05695   0.17388   0.01700;
        6   3   4   0.06701   0.17103   0.01730;
        7   4   5   0.01335   0.04211   0.00640;
        8   4   7   0         0.20912   0;
        9   4   9   0         0.55618   0;
        10  5   6   0         0.25202   0;
        11  6   11  0.09498   0.19890   0;
        12  6   12  0.12291   0.25581   0;
        13  6   13  0.06615   0.13027   0;
        14  7   8   0         0.17615   0;
        15  7   9   0         0.11001   0;
        16  9   10  0.03181   0.08450   0;
        17  9   14  0.12711   0.27038   0;
        18  10  11  0.08205   0.19207   0;
        19  12  13  0.22092   0.19988   0;
        20  13  14  0.17093   0.34802   0;
    ];
    
    % Transformer tap data from Table A.3
    tap_data = [
        4   7   0.978;
        4   9   0.969;
        5   6   0.932;
    ];
    
    % Shunt capacitor data from Table A.5
    shunt_data = [
        9   0.19;
    ];
    
    % Base values
    Sbase = 100; % MVA
    Vbase = 1;   % pu (assuming all voltages are in pu)
    
    % Extract bus information
    bus_number = bus_data(:,1);
    V_mag = bus_data(:,2);
    V_angle = bus_data(:,3);
    P_gen = bus_data(:,4) / Sbase;
    Q_gen = bus_data(:,5) / Sbase;
    P_load = bus_data(:,6) / Sbase;
    Q_load = bus_data(:,7) / Sbase;
    Q_min = bus_data(:,8) / Sbase;
    Q_max = bus_data(:,9) / Sbase;
    
    % Number of buses
    nbus = length(bus_number);
    
    % Initialize Ybus
    Ybus = zeros(nbus, nbus);
    
    % Build Ybus matrix
    for k = 1:size(line_data,1)
        from_bus = line_data(k,2);
        to_bus = line_data(k,3);
        R = line_data(k,4);
        X = line_data(k,5);
        B = line_data(k,6);
        
        % Check if this line is a transformer
        is_transformer = 0;
        for m = 1:size(tap_data,1)
            if (from_bus == tap_data(m,1) && to_bus == tap_data(m,2)) || ...
               (from_bus == tap_data(m,2) && to_bus == tap_data(m,1))
                tap = tap_data(m,3);
                is_transformer = 1;
                break;
            end
        end
        
        if is_transformer
            % Transformer admittance
            y = 1 / (R + 1j*X);
            Ybus(from_bus, to_bus) = Ybus(from_bus, to_bus) - y/tap;
            Ybus(to_bus, from_bus) = Ybus(to_bus, from_bus) - y/tap;
            Ybus(from_bus, from_bus) = Ybus(from_bus, from_bus) + y/(tap^2);
            Ybus(to_bus, to_bus) = Ybus(to_bus, to_bus) + y;
        else
            % Regular line admittance
            y = 1 / (R + 1j*X);
            Ybus(from_bus, to_bus) = Ybus(from_bus, to_bus) - y;
            Ybus(to_bus, from_bus) = Ybus(to_bus, from_bus) - y;
            Ybus(from_bus, from_bus) = Ybus(from_bus, from_bus) + y + 1j*B/2;
            Ybus(to_bus, to_bus) = Ybus(to_bus, to_bus) + y + 1j*B/2;
        end
    end
    
    % Add shunt capacitors
    for k = 1:size(shunt_data,1)
        bus = shunt_data(k,1);
        B_shunt = shunt_data(k,2);
        Ybus(bus, bus) = Ybus(bus, bus) + 1j*B_shunt;
    end
    
    % Determine bus types
    % 1 = Slack bus (bus 1)
    % 2 = PV bus (buses 2, 3)
    % 3 = PQ bus (all others)
    bus_type = ones(nbus,1) * 3; % Initialize all as PQ
    bus_type(1) = 1; % Slack bus
    bus_type(2) = 2; % PV bus
    bus_type(3) = 2; % PV bus
    
    % Initialize voltages
    V = V_mag .* exp(1j * V_angle * pi/180);
    
    % Calculate net injected power
    P_spec = P_gen - P_load;
    Q_spec = Q_gen - Q_load;
    
    % Gauss-Seidel parameters
    max_iter = 100;
    tol = 1e-6;
    alpha = 1.6; % Acceleration factor (1.0 for no acceleration)
    
    % Start Gauss-Seidel iterations
    for iter = 1:max_iter
        V_prev = V;
        
        % Update voltages bus by bus
        for i = 1:nbus
            if bus_type(i) == 1
                % Slack bus - voltage stays constant
                continue;
            end
            
            % Calculate current injection
            S = P_spec(i) + 1j*Q_spec(i);
            I = conj(S / V(i));
            
            % Calculate new voltage estimate
            sum_YV = 0;
            for k = 1:nbus
                if k ~= i
                    sum_YV = sum_YV + Ybus(i,k) * V(k);
                end
            end
            
            V_new = (I - sum_YV) / Ybus(i,i);
            
            if bus_type(i) == 2
                % PV bus - keep magnitude constant, update angle
                V_new = V_mag(i) * exp(1j * angle(V_new));
                
                % Calculate reactive power injection
                Q_calc = -imag(conj(V_new) * (Ybus(i,:) * V));
                
                % Check reactive power limits
                if Q_calc < Q_min(i)
                    Q_calc = Q_min(i);
                    bus_type(i) = 3; % Switch to PQ bus if limit violated
                elseif Q_calc > Q_max(i)
                    Q_calc = Q_max(i);
                    bus_type(i) = 3; % Switch to PQ bus if limit violated
                end
                
                % Update Q_spec with calculated value
                Q_spec(i) = Q_calc;
            end
            
            % Apply acceleration factor
            V(i) = V_prev(i) + alpha * (V_new - V_prev(i));
        end
        
        % Check for convergence
        max_diff = max(abs(V - V_prev));
        if max_diff < tol
            fprintf('Converged in %d iterations.\n', iter);
            break;
        end
        
        % Display progress every 50 iterations
        if mod(iter, 50) == 0
            fprintf('Iteration %d: Max voltage change = %.6f pu\n', iter, max_diff);
        end
    end
    
    if iter == max_iter
        fprintf('Did not converge in %d iterations. Max voltage change = %.6f pu\n', max_iter, max_diff);
    end
    
    % Display results
    fprintf('\nBus Voltage Results:\n');
    fprintf('Bus #   Voltage (pu)   Angle (deg)   P Gen (MW)   Q Gen (MVAR)   P Load (MW)   Q Load (MVAR)\n');
    for i = 1:nbus
        fprintf('%4d %12.4f %12.4f %12.2f %12.2f %12.2f %12.2f\n', ...
                bus_number(i), abs(V(i)), angle(V(i))*180/pi, ...
                P_gen(i)*Sbase, Q_gen(i)*Sbase, P_load(i)*Sbase, Q_load(i)*Sbase);
    end
    
    % Calculate line flows
    fprintf('\nLine Power Flows:\n');
    fprintf('From Bus   To Bus   P Flow (MW)   Q Flow (MVAR)\n');
    for k = 1:size(line_data,1)
        from_bus = line_data(k,2);
        to_bus = line_data(k,3);
        R = line_data(k,4);
        X = line_data(k,5);
        B = line_data(k,6);
        
        % Check if this line is a transformer
        is_transformer = 0;
        tap = 1;
        for m = 1:size(tap_data,1)
            if (from_bus == tap_data(m,1) && to_bus == tap_data(m,2)) || ...
               (from_bus == tap_data(m,2) && to_bus == tap_data(m,1))
                tap = tap_data(m,3);
                is_transformer = 1;
                break;
            end
        end
        
        if is_transformer
            if from_bus == tap_data(m,1) % From is primary side
                I_from = (V(from_bus)/tap - V(to_bus)) / (R + 1j*X);
                S_from = V(from_bus)/tap * conj(I_from);
            else % From is secondary side
                I_from = (V(from_bus) - V(to_bus)/tap) / (R + 1j*X);
                S_from = V(from_bus) * conj(I_from);
            end
            I_to = -I_from;
            S_to = V(to_bus) * conj(I_to);
        else
            I_from = (V(from_bus) - V(to_bus)) / (R + 1j*X) + V(from_bus) * (1j*B/2);
            S_from = V(from_bus) * conj(I_from);
            I_to = (V(to_bus) - V(from_bus)) / (R + 1j*X) + V(to_bus) * (1j*B/2);
            S_to = V(to_bus) * conj(I_to);
        end
        
        fprintf('%4d %9d %12.2f %12.2f\n', from_bus, to_bus, real(S_from)*Sbase, imag(S_from)*Sbase);
    end
end

