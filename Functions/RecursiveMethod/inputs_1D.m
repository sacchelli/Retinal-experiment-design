function [controls, counter, cStar, kMax, inputParameters] = inputs_1D(m, dn, N, Delta, Tf, retinalWidth,...
                                addRandomInputs, nbrRandomInputs, seed )
% function [controls, counter, cStar, kMax, inputParameters] = inputs_1D(m, dn, N, Delta, Tf, retinalWidth,...
%                                addRandomInputs, nbrRandomInputs, seed )



fprintf(['\n','Input collection... '])

tStepStart = tic;

cStar = retinalWidth/(100*Delta);
% cStar = retinalWidth/(Tf);

cMax = 20*cStar;   % max wave velocity
cMin = -20*cStar; % min wave velocity
% cMin = 0; % min wave velocity

kMax = 2*pi/retinalWidth*(dn-1)/2;   % max wave number
kMin = 2*pi/(retinalWidth);   % min wave number

cN = 51;
kN = 30;

[controls1DWave, inputParameters1Dwave] = generate1DPlaneWaveControls(dn, N, Delta, Tf, retinalWidth, cMin, cMax, cN, kMin, kMax, kN);

if addRandomInputs % produce random inputs if asked
    [controlsRandom, inputParametersRandom] = generate1DRandomControls(dn, N, nbrRandomInputs, seed);
end


% Concatenate inputs and constant input.
if addRandomInputs 
    counter = length(controls1DWave) + length(controlsRandom) + 1;
    inputParameters = [inputParameters1Dwave, inputParametersRandom, [NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = [controls1DWave; controlsRandom];
else
    counter = length(controls1DWave) + 1;
    inputParameters = [inputParameters1Dwave,  [NaN;0]];
    controls = cell(counter,1);
    controls(1:end-1) = controls1DWave;
end

%

controls{counter} = ones(m, N); % Constant input at the end
K = counter;
tStepEnd = toc(tStepStart);

fprintf(['done (',num2str(tStepEnd,3) ,'s) \n'])

fprintf(['There are ', num2str(K),' inputs in this experiment'])

end