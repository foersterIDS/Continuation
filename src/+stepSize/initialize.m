%% path continuation - stepSize.initialize
%
%   Institute of Dynamics and Vibration Research
%   Leibniz University Hannover
%   18.01.2022 - Tido Kubatschek
%
function [ds0] = initialize(oih,var0,lStart,lEnd,ds0)
    %% check chosen stepsize options
    %
    if oih.opt.checkStepSizeOptions
        %
        if ds0 < oih.opt.dsMin % ds0 must not be lower than dsMin
            ds0 = oih.opt.dsMin;
            aux.printLine(oih,'--> ds0 is too small! It must not be smaller than dsMin. ds0 has been adapted.\n');
        end
        %
        if ds0 > oih.opt.dsMax % ds0 must not be larger than dsMax
            ds0 = oih.opt.dsMax;
            aux.printLine(oih,'--> ds0 is too large! It must not be greater than dsMax. ds0 has been adapted.\n');
        end
        %
        if oih.opt.dsMin < 10*oih.opt.solverTol % dsMin cannot be lower than tolerance of solver
            oih.opt.dsMin = 10*oih.opt.solverTol;
            aux.printLine(oih,'--> dsMin has to be at least 10*solverTol = %.2e. dsMin has been adapted.\n',oih.opt.dsMin);
        end
        %
        if ds0 > (sqrt(sum(var0.^2) + (lEnd-lStart)^2)/10) % dsMax should not be larger than mag of points
            % get order of magnitude
            nMag = floor(log10(sqrt(sum(var0.^2) + (lEnd-lStart)^2)));
            % adapt ds0
            ds0 = 10^(nMag-1);
            % order of magnitude must be at least 1 lower
            aux.printLine(oih,'--> ds0 is too large! It must not be greater than 1.00e%i. ds0 has been adapted.\n',nMag-1);
        end
        %
    end
    %
    %% Tell StepsizeOptions which values must be calculated
    %
    if oih.opt.stepSizeControl.multiplicative || oih.opt.stepSizeControl.multiplicativeAlt
        %% check weigths
        %
        weights = oih.opt.weightsMultiplicative;
        stepsizeOptions = struct('speedOfContinuation',weights(3) ~= 0,...
                                  'predictor', weights(5) ~= 0,...
                                  'rateOfContraction',weights(4) ~= 0,...
                                  'iterations', true);
    elseif oih.opt.stepSizeControl.error || oih.opt.stepSizeControl.errorAlt
        %% check weigths
        %
        weights = oih.opt.weightsError;
        stepsizeOptions = struct('speedOfContinuation',weights(3) ~= 0,...
                                  'predictor', weights(5) ~= 0,...
                                  'rateOfContraction',weights(4) ~= 0,...
                                  'iterations', true);
    elseif oih.opt.stepSizeControl.contraction
        stepsizeOptions = struct('speedOfContinuation',false,...
                                  'predictor',false,...
                                  'rateOfContraction',true,...
                                  'iterations', true);
    elseif oih.opt.stepSizeControl.fix
        stepsizeOptions = struct('speedOfContinuation',false,...
                                  'predictor',false,...
                                  'rateOfContraction',false,...
                                  'iterations', true);
        oih.opt.dsMax = ds0;
        oih.opt.dsMin = ds0;
    else
        stepsizeOptions = struct('speedOfContinuation',false,...
                                  'predictor',false,...
                                  'rateOfContraction',false,...
                                  'iterations', true);
    end
    %
    %% set step size options in oih
    %
    oih.stepsizeOptions = stepsizeOptions;
    %
end