function [] = quailsMinimizer(fileID)

    % default file
    if nargin < 1
        fileID = '';
    end

    addpath('constraints');
    addpath('meritfunctions');
    addpath('guesses');

    % display info
    disp('SHOTGUN MINIMIZATION for PERTURBATIVE CHROMATICITY CORRECTION')
    disp('Searches parameter space using shotgun approach.')

    % minimisation parameters and options
    [lb, ub, maxStepShots, maxShots, funEvals, funTol, xTol] = matchingParams();
    options = optimoptions(@fmincon, 'Algorithm', 'sqp', 'TolX', xTol, 'TolFun', funTol, 'MaxFunEvals', funEvals);

    % display parameters
    disp(['--- Tolerance: ' num2str(funTol) ' (func.), ' num2str(xTol) ' (params), Evals: ' num2str(funEvals) ' ---']);
    disp(['--- Max attempts: ' num2str(maxShots) ' (total), ' num2str(maxStepShots) ' (same step) ---']);
    disp(['--- Upper limits: ' mat2str(ub) ' ---']);
    disp(['--- Lower limits: ' mat2str(lb) ' ---']);

    % improvement variables (saves the best solution)
    minfval = Inf;
    bestx = [];
    stepshot = 1;

    guessFlag = false;
    evalc(['guess = initialGuess_' fileID '()']);
    for shot = 1:maxShots

        % make random guess (shotgun approach)

        if(numel(guess) == numel(lb))
            x0 = guess;
            guessFlag = true;
        else
            x0 = (ub-lb).*rand(1,length(ub))+lb;
        end

        % silent evaluation of fmincon
        evalc(['[x, fval,exitflag] = fmincon(@meritfunction_' fileID ', x0, [],[],[],[], lb, ub, @constraints_' fileID ', options);']);
        %disp(['Trying: ' num2str(x0)]);

        % write out if better results
        evalc(['[~, consts] = constraints_' fileID '(x)']);
        constsum = sum(consts.*consts);
        matched = (constsum < 1e-4);

        if(fval < 0.99*minfval && fval >= 0 && matched)
            minfval = fval;
            bestx = x;
            disp(['> ' num2str(stepshot) ' attempts made.']);
            disp(['Min. value: ' num2str(minfval,'%2.2e') ' -- ' mat2str(bestx,3) ' (match: ' num2str(constsum) ')']);
            saveSolution(bestx, fileID);
            stepshot = 0;
        end

        % break if too many shots on same step
        if stepshot >= maxStepShots
            break;
        end

        stepshot = stepshot + 1;

        % exit if guess
        if(numel(guess) == numel(lb))
            break;
        end

    end
    
    disp(['> ' num2str(stepshot) ' unsuccessful attempts made.']);
    disp(' ');
    disp(['BEST SOLUTION: ' num2str(minfval,'%2.2e') ' : ' mat2str(bestx,4)])

    % end if guess (no point doing it again)
    if guessFlag
        return;
    end
end
