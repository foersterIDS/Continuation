%% path continuation - stepSize.initialize
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.01.2022 - Tido Kubatschek
%
function [Opt,ds0,StepsizeOptions] = initialize(Opt,var0,lStart,lEnd,ds0)
    %% check chosen stepsize options
    %
    if Opt.checkStepSizeOptions
        %
        if ds0 < Opt.dsMin % ds0 must not be lower than dsMin
            ds0 = Opt.dsMin;
            aux.printLine(Opt,'--> ds0 is too small! It must not be smaller than dsMin. ds0 has been adapted.\n');
        end
        %
        if ds0 > Opt.dsMax % ds0 must not be larger than dsMax
            ds0 = Opt.dsMax;
            aux.printLine(Opt,'--> ds0 is too large! It must not be greater than dsMax. ds0 has been adapted.\n');
        end
        %
        if Opt.dsMin < 10*Opt.solverTol % dsMin cannot be lower than tolerance of solver
            Opt.dsMin = 10*Opt.solverTol;
            aux.printLine(Opt,'--> dsMin has to be at least 10*solverTol = %.2e. dsMin has been adapted.\n',Opt.dsMin);
        end
        %
        if ds0 > (sqrt(sum(var0.^2) + (lEnd-lStart)^2)/10) % dsMax should not be larger than mag of points
            % get order of magnitude
            nMag = floor(log10(sqrt(sum(var0.^2) + (lEnd-lStart)^2)));
            % adapt ds0
            ds0 = 10^(nMag-1);
            % order of magnitude must be at least 1 lower
            aux.printLine(Opt,'--> ds0 is too large! It must not be greater than 1.00e%i. ds0 has been adapted.\n',nMag-1);
        end
        %
    end
    %
    %% Tell StepsizeOptions which values must be calculated
    
    if Opt.stepSizeControl.multiplicative || Opt.stepSizeControl.multiplicativeAlt
        %% check weigths
        %
        weights = Opt.weightsMultiplicative;
        StepsizeOptions = struct('speedOfContinuation',weights(3) ~= 0,...
                                  'predictor', weights(5) ~= 0,...
                                  'rateOfContraction',weights(4) ~= 0,...
                                  'iterations', true);
    elseif Opt.stepSizeControl.error || Opt.stepSizeControl.errorAlt
        %% check weigths
        %
        weights = Opt.weightsError;
        StepsizeOptions = struct('speedOfContinuation',weights(3) ~= 0,...
                                  'predictor', weights(5) ~= 0,...
                                  'rateOfContraction',weights(4) ~= 0,...
                                  'iterations', true);
    elseif Opt.stepSizeControl.contraction
        StepsizeOptions = struct('speedOfContinuation',false,...
                                  'predictor',false,...
                                  'rateOfContraction',true,...
                                  'iterations', true);
    elseif Opt.stepSizeControl.fix
        StepsizeOptions = struct('speedOfContinuation',false,...
                                  'predictor',false,...
                                  'rateOfContraction',false,...
                                  'iterations', true);
        Opt.dsMax = ds0;
        Opt.dsMin = ds0;
    else
        StepsizeOptions = struct('speedOfContinuation',false,...
                                  'predictor',false,...
                                  'rateOfContraction',false,...
                                  'iterations', true);
    end
    %
end