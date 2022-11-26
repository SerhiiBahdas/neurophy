function nY = impEuler(odefun, nY0, tTime, tStep, varargin)
    % IMPEULER Solve system of ordinary differential equations using
    %   implicit Euler method. 
    %
    %   nY = impEuler(odefun, nY0, tTime, tStep, varargin)
    %
    %   INPUT ===========================================================
    %
    %   odefun(function)
    %   Specifies a system of differential eqautions. 
    %   Example: odefun = @(t,y) [-100*y(1); -200*y(2)]
    %
    %   nY0 (numeric array)
    %   Initial conditions. 
    %   Example: nY0 = [1,1]
    %
    %   tTime (numeric array)
    %   Interval of integration. 
    %   Example: tTime = [0,0.2];
    %
    %   tStep (double)
    %   Integration step. 
    %   Example: tStep = 1e-3; 
    %
    %   bPlot (boolean)
    %   Visualize the solution.
    %   Example: 1
    %
    %   OUTPUT ==========================================================
    %
    %   nY (numeric array)
    %   Solution of the system of ODE. 
    %   
    %   EXAMPLE =========================================================
    %
    %   odefun = @(t,y) [-100*y(1); -200*y(2)];
    %   nY0    = [1,1]; 
    %   tTime  = [0,0.2];
    %   tStep  = 1e-3; 
    %   bPlot  = 1; 
    %   nY = impEuler(odefun, nY0, tTime, tStep, 'bPlot', 1);
    %
    %   AUTHOR ==========================================================
    %
    %   S.Bahdasariants, NEL, WVU, https://github.com/SerhiiBahdas
    %
    %   =================================================================
    
    %% ASSIGN INPUT VARIABLES TO A STRUCTURE
    
    % Check number of function input arguments 
    if numel(varargin)/2 ~= 0 
        
        % Create a structure where the parameters will be stored
        in = struct(); 
        
        % For all additional arguments 
        for iInput = 1:numel(varargin)/2
            
            % Determine indexes of the variable name
            %   and value
            idxVarName  = (iInput-1)*2 + 1; 
            idxVarValue = (iInput-1)*2 + 2; 
            
            % Assign argument's name and value to a struct
            in.(varargin{idxVarName}) = varargin{idxVarValue};
            
        end % iInput
       
    % Otherwise, assign default parameters    
    else
        
        % Opt not to visualize solutions
        in.bPlot = 0; 
        
    end % if
    
    %% CHECK INPUTS VALIDITY
    
    % Determine integration time interval 
    tSpan = tTime(2) - tTime(1); 

    % If the time vector was specified incorrectly
    if tSpan <=0, error('Specify valid integration interval.'); end

    % Determine number of itterations 
    nItter = round(tSpan/tStep); 
    
    % If the number of itterations is less than 1, throw an error
    if floor(tSpan/tStep)<1, error('Choose smaller integration step.'); end
    
    % Preallocate memory for solution
    nY = zeros(nItter, numel(nY0)); 
    
    % Reshape initial consitions 
    nY0 = reshape(nY0,[],2)'; 
    
    %% SOLVE SYSTEM OF ODE
    
    % Find solution of the ODE
    for iItter = 1:nItter

        % Time of the i+1 sample
        nTime = tTime(1) + iItter*tStep; 

        % Implicit Euler equation
        imEuler = @(Y) Y - nY0 - tStep*odefun(nTime,Y);

        % Use Netwon's method to solve nonlin. eq. 
        nY(iItter,:) = quaziNewton(imEuler, nY0); 

        % Update the previous value of the function 
        nY0 = nY(iItter,:); 
        
        % Reshape the solution
        nY0 = reshape(nY0,[],2)'; 

    end % iItter

    %% VISUALIZE SOLUTION
    
    % If the results need to be visualized 
    if in.bPlot == 1
        
        figure(); 
        
        % Create a time vector 
        nTime = linspace(tTime(1), tTime(2), nItter); 
        
        % Plot the solution(s)
        plot(nTime, nY, "LineWidth", 3);
        
        % Add labels
        xlabel('t'); 
        ylabel('y'); 
        
    end % if
    
end % impEuler

function nY = quaziNewton(fun, nY0, varargin)
    % QUAZINEWTON Solves the system of nonlinear equations
    %   using quzi-Newton's method.
    %
    %   INPUT ============================================
    %
    %   fun (function)
    %   System of nonlinear equations to solve.
    %   Example: [x(1) + x(2) - 3; x(1).^2 + x(2).^2 - 9]
    %
    %   nY0 (numeric array)
    %   Estimation of the solution. 
    %   Example [0.01, 2.98]
    %
    %   numStep (double)
    %   Maximum number of itterations. 
    %   Example: 100
    %
    %   maxError (double)
    %   Maximum acceptable error.
    %   Example: 1e-6; 
    %
    %   OUTPUT ===========================================
    %
    %   nY (numeric array)
    %   Solution of the system of nonlinear equations. 
    %
    %   EXAMPLE ==========================================
    %
    %   fun = @(x) [x(1) + x(2) - 3; x(1).^2 + x(2).^2 - 9];
    %   nY0 = [0.01, 2.98];
    %   nY = quaziNewton(fun, nY0, 'numStep', 100, 'maxError', 1e-6);
    %
    %   AUTHOR ===========================================
    %
    %   S.Bahdasariants, NEL, WVU, serhiibahdas@gmail.com
    %
    %   ==================================================
    
    % Check number of function input arguments 
    if numel(varargin)/2 ~= 0 
        
        % Create a structure where the parameters will be stored
        in = struct(); 
        
        % For all additional arguments 
        for iInput = 1:numel(varargin)/2
            
            % Determine indexes of the variable name
            %   and value
            idxVarName  = (iInput-1)*2 + 1; 
            idxVarValue = iInput*2; 
            
            % Assign argument's name and value to a struct
            in.(varargin{idxVarName}) = varargin{idxVarValue};
            
        end % iInput
       
    % Otherwise, assign default parameters    
    else
        
        % Assign the maximum number of steps 
        in.numStep = 100; 
        
        % Define maximum acceptable error
        in.maxError = 1e-6; 
        
    end % if
    
    % Reshape estimation of the solution
    nY0 = reshape(nY0,[],2)'; 
    
    % Find suboptimal step used to estimate a derivative of
    %   the nonlinear equations. 
    
    nStep = sqrt(eps)*abs(nY0(1)); 
    
    %% ESTIMATE JACOBIAN USING FINITE DIFFERENCES METHOD
    
    % Preallocate memory for Jacobian matrix 
    J0 = zeros(numel(nY0), numel(nY0)); 
    
    % Compute partial derivatives with respect to the
    %   i-th variable to compose Jacobian column by 
    %   column. 
    
    % For all variables in the equations
    for iVar = 1:numel(nY0)

        % Increase currently analyzed variable by 
        %   small value to later compute a derivative
        
        nY_eps = nY0; 
        nY_eps(iVar) = nY_eps(iVar) + nStep; 
        
        % Compute i-th column of the Jacobian using 
        %   finite differences method
        
        J0(:,iVar) = (fun(nY_eps) - fun(nY0))/nStep; 

    end % iVar
        
    % Find the error of the first estimate of the 
    %   solution
    nDeltaY0 = pinv(J0)*(-fun(nY0)); 
    
    %% QUAZI-NEWTON'S METHOD
    
    % Create a function calculating Eucledian norm
    euclid = @(y) sqrt(sum(y(:).*y(:)));
    
    % Run quazi-Newton method to find the roots of the
    %   nonlinear eqaution
    
    % For the number of itterations specified
    for iItter = 1:in.numStep
        
        % Compute solutions for the i-th step
        nY1 = nY0 + nDeltaY0; 
        
        % Plug in the solutions on the i-th 
        %   step into the system of equations
        nF1 = fun(nY1); 
        
        % Compute Jacobin using Broyden's method
        J1 = J0 + nF1*transpose(nDeltaY0)/...
            (transpose(nDeltaY0)*nDeltaY0);
        
        % Compute the error of the estimate 
        nDeltaY1 = pinv(J1)*(-nF1);
        
        % Check if the solution satisfies criteria
        if euclid(nDeltaY1)/euclid(nY1) <= in.maxError
            break;
        end % if
        
        % Reassign the variables for the next itteration
        J0 = J1; nY0 = nY1; nDeltaY0 = nDeltaY1; 

    end % iItter
    
    if euclid(nDeltaY1)/euclid(nY1) > in.maxError 
    	error('The solution was not found. Relax maximum acceptable error or increase number of itterations.'); 
    end % if
    
    % Save the solution
    nY = nY1; 
    
end % quaziNewton
