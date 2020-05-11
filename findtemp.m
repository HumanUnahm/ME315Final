% Manually clear to restart program
if(~exist('workspaceCleared'))
    clc, clear
    workspaceCleared=1;
    prevNodes = 0;
    prevdt = 0;
    prevtf = 0;
    runSim = 0;
end

% Which Problem do you want solved? (Set all other options to false)
Req4_Nodes101 = false;
Req4_Nodes201 = true;
Req5 = false;
Req6 = false;
Req7 = false;

% Update Simulation Settings
if(Req4_Nodes101||Req5)
    Nodes = 101; % Number of Nodes
    dt = 0.025; % Time Step [s]
    tf = 60*5; % Simulation Length [s]
end

if(Req4_Nodes201||Req5&&~(Req4_Nodes101||Req5))
    Nodes = 202; % Number of Nodes
    dt = 0.010; % Time Step [s]
    tf = 60*5; % Simulation Length [s]
end
if(Req5)
    Nodes = 101; % Number of Nodes
    dt = 0.025; % Time Step [s]
    tf = 60*60; % Simulation Length [s]
end

if(Req6)
    Nodes = 101; % Number of Nodes
    dt = 0.025; % Time Step [s]
    tf = 60*60; % 1 hour length
end
if(Req7)
    Nodes = 101; % Number of Nodes
    dt = 0.025; % Time Step [s]
    tf = 7000; % 7000 seconds (Refenced gifted solutions)
end
%% Stops simulation from being repeated using same sim parameters
if(prevNodes ~= Nodes || prevdt ~= dt || prevtf ~= tf)
    prevNodes = Nodes;
    prevdt = dt;
    prevtf = tf;
    runSim = true;
else
    runSim = false;
end



%% Simulation & System Parameters
if(runSim)
            % All units must be base SI [kg,W,J,m]
    %Analytical Parameters
    eigenValueLimit = 50;
    err = 0.001; % Decrease if getting repeat eigen values.
    % Room Properties
    h = 150; % Convection Heat Transfer Coeffecient [W/m^2/K]
    T_amb = 24; % Room Ambient Temperature [Celsius]
    % Sphere Properties (Steel)
    Cp = 510; % Heat Capacity [j/kg/K]
    rho = 7810; % Density [kg/m^3]
    k = 11;  % Conduction Heat Transfer Coeffecient [W/m/K]
    r = 0.1;  % Sphere radius [m]
    T_i = 350; % Sphere Initial TemperatureTemperature [Celsius]
    
    %% Initialization Part 1 [Intermediate Variables and Allocation]
    A = zeros(Nodes,Nodes); % Allocate Coeffecient Matrix
    b = zeros(Nodes,1); % Allocate Constants Vector
    r_pos = linspace(0,r,Nodes); % position(node) Vector
    dr = r_pos(2); % radial step (Assumes r_vec(1) = 0)
    alpha = k/rho/Cp; % Heat Transport Coeffecient
    t_vec = 0:dt:tf; % Time Vector
    I = eye(Nodes); % Identity Matrix for matrix [A].
    NumOfStates = length(t_vec); % Number of States
    T = zeros(Nodes,NumOfStates); % Allocate Temperature matrix; Columns are states, rows are nodes.
    T_ana = zeros(Nodes,1); % Allocate Analytical Temperature matrix; Columns are states, rows are nodes.
    lambda = zeros(eigenValueLimit,1); % Allocate Eigenvalues
    
    %% Initialization Part 2 [Coeffecient Matrix, Initial State, Constant Matrix]
    T(:,1) = T_i; % Set all nodes of first state to initial temperature

    % A5 appendix equation for first row (T1 and T2 nodes only)
    A(1,1) = -6*alpha/dr^2; %T_1 coeffecient
    A(1,2) = 6*alpha/dr^2; %T_2 coeffecient
    b(1) = 0; % A5 has no constant

    % Eq6 for interior nodes (start after first node, end before last node)
    for i = 2:Nodes-1
        A(i,i-1) = alpha*(r_pos(i)*dr-dr^2)/(r_pos(i)*dr^3);
        A(i,i) = -2*alpha/dr^2;
        A(i,i+1) = alpha*(dr^2+r_pos(i)*dr)/(r_pos(i)*dr^3);
        b(i) = 0;% Interior nodes have no constant
    end
    % A6 appendix equation for last row(T_N and T_N-1 nodes).
    A(Nodes,Nodes-1) = alpha*2/dr^2; % 2nd before last node
    A(Nodes,Nodes) = alpha*(-2*h/r/k-2*h/k/dr-2/dr^2); % Last Node
    b(Nodes) = alpha*(2*h/k/r+2*h/k/dr)*T_amb; % Last Constant

    %% Simulate
    % Run Inverse Calculation only once.
    M1 = (I-dt/2*A)\I;
    M2 = (I+dt/2*A);

    for j = 1:NumOfStates-1 % Use j notation now for system states
        T(:,j+1) = M1*M2*T(:,j)+dt*M1*b;% Eq 7
    end
end



%% Analytical Solution Stuff

% function handle for numerical solve function fzero.
f = @(lambda) lambda*cos(lambda)+(h*r/k-1)*sin(lambda);
% obtain first eigein value (after zero)
z = fzero(f,1);
lambda(1) = z;

i=1; % Eigen Index
x=1; % Domain variable
while ( i < eigenValueLimit)
    x = x + 1; % Step through x-domain by 1.
    z = fzero(f,x); % Find next zero
    if (abs(lambda(i)-z) < err) % If zero is the same (within bounds) start over
        continue;
    end
    % Zero found sucesfully, add it to the list of eigen values.
    i = i + 1;
    lambda(i) = z; 
end


for i=1:Nodes
    r_dist = r_pos(i); % Grab Position of node
    t_state = 60*5; % 5mins in seconds
    total = 0; % Reset Sum
    for j=1:eigenValueLimit
        lamb = lambda(j);
        numerator = (sin(lamb)-lamb*cos(lamb))*sin(lamb*r_dist/r)...
            *exp(-lamb^2*alpha*t_state/r^2);
        deonominator = (2*lamb-sin(2*lamb))*lamb*r_dist/r;
        total = total + numerator/deonominator;
    end
    T_ana(i)=4*total*(T_i-T_amb)+T_amb;
end



%% Questions that need answering
% Analytical Simulation Comparison (Low Quality)
if(Req4_Nodes101)
    % Eigen Value vs Index N
    figure; hold on; grid on;
    plot(1:1:eigenValueLimit, lambda)
    title('Eigen Value vs Index n');
    ylabel('EigenValue'); 
    xlabel('Index n');
    
    % Calculate RMSE And Percent Error
    stepspermin = 60/dt; % How many steps per minute?
    time_index = 5*stepspermin; % Index for time @ 5minutes
    TempState = T(:,time_index);% Temperature State of System @ 5minutes
    RMSE = sqrt(sum((T_ana-TempState).^2)/Nodes);
    PercentError = abs(T_ana-TempState)./T_ana*100;
    figure; hold on; grid on;
    plot(r_pos*1e2, PercentError);
    MaxPercentError101Nodes = max(PercentError);
    xlabel('r[cm]'); 
    ylabel('Absolute Percent Error[%]');
    title('101 nodes; @ 5 mins and dt=0.025secs');
end
% Analyatical Simulation Comparison (High Quality)
if(Req4_Nodes201)
    % Eigen Value vs Index N
    figure; hold on; grid on;
    plot(1:1:eigenValueLimit, lambda)
    title('Eigen Value vs Index n');
    ylabel('EigenValue'); 
    xlabel('Index n');
    
    % Calculate RMSE And Percent Error
    stepspermin = 60/dt; % How many steps per minute?
    time_index = 5*stepspermin; % Index for time @ 5minutes
    TempState = T(:,time_index);% Temperature State of System @ 5minutes
    RMSE = sqrt(sum((T_ana-TempState).^2)/Nodes);
    PercentError = abs(T_ana-TempState)./T_ana*100;
    figure; hold on; grid on;   
    plot(r_pos*1e2, PercentError);
    MaxPercentError201Nodes = max(PercentError);
    xlabel('r[cm]'); 
    ylabel('Absolute Percent Error[%]');
    title('201 nodes; @ 5 minutes and dt=0.01secs');
end
% Simulation Isotime graphs for 60secs and 1 hour
if(Req5)
    % Data has been generated. All that's left to do is graph.
    figure; hold on; grid on;
    stepsPerSec = 1/dt;
    timeplot = [1,1,2,5,10,20,30,45,60]*stepsPerSec;
    timeplot(1)=1;
    plot(r_pos*1e2, T(:,timeplot))
    xlabel('x - [cm]'); 
    ylabel('T - [C]');
    legend('0sec','1sec','2sec','5sec','10sec','30sec','45sec','60sec')
    title('101 nodes, timestep 25ms, 60 seconds simulation')

    figure; hold on; grid on;
    stepspermin = 60/dt;
    timeplot = [1,1,2,5,10,20,30,45,60]*stepspermin;
    timeplot(1)=1;
    plot(r_pos*1e2, T(:,timeplot))
    xlabel('r[cm]'); 
    ylabel('T[C]');
    legend('0min','1min','2min','5min','10min','20min','30min','45min','60min');
    title('101 nodes, timestep 25ms, 1 hour simulation')
end
%% Stuff thats more straight forward and don't need allocation/initialization.
%Volume Averaged Temperature at 1 hour range
if(Req6)
    T_numerical = zeros(NumOfStates,1);
    for j=1:NumOfStates % Iterate through states and apply volume avergte temperature formula
        T_numerical(j) = (3/r^3*trapz(r_pos, T(:,j)'.*(r_pos.^2))); % Integration done using trapzoidal summation
    end
    % Lumped Capacitence
    V = 4/3*pi*r^3;
    A = 4*pi*r^2;
    CharacteristicLength = V/A;
    Biot = h*CharacteristicLength/k;
    T_capacitence = zeros(NumOfStates,1);
    for j=1:NumOfStates
        FourierNumber = alpha*t_vec(j)/CharacteristicLength^2; % Calc Fouier Number
        T_capacitence(j) = T_amb+(T_i-T_amb)*exp(-Biot*FourierNumber); % Apply Lumped Capacitence Formula
    end
    figure; hold on; grid on;
    plot(t_vec/60, T_numerical)
    plot(t_vec/60, T_capacitence)
    title('Numerical Average Temperature VS Lumped Capacitance Model');
    legend('Volume Averaged (Numerical)', 'Lumped Capacitance')
    xlabel('Time[mins]');
    ylabel('T[C]');
end

%% Volume Averaged through 7000secs (Copy Pasta from above)
if(Req7)
    T_numerical = zeros(NumOfStates,1);
    for j=1:NumOfStates % Iterate through states and apply volume averge temperature formula
        % Calculate Temperature, state by state, be giving trapz temperature
        % states
        T_state = T(:,j)';
        T_numerical(j) = (3/r^3*trapz(r_pos, T_state.*(r_pos.^2))); % Integration done using trapzoidal summation
    end
 
    figure; hold on; grid on;
    plot(t_vec/60, T_numerical)
    title('Numerical Average Temperature');
    legend('Volume Averaged (Numerical)')
    xlabel('Time[mins]');
    ylabel('T[C]');
end
% Out of time :(