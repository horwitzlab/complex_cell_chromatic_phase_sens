% PhaseNTsim.m
% Simulation of an experiment to probe the interactions between chromatic
% and achromatic contributions to the responses of color-complex cells.
%
% Simulated neurons
% All neurons in this simulation have noiseless, real-valued responses.
% The simulated color-complex cell receives input from two pairs of two 
% linear neurons (2 color directions x cos and sin phases). Each linear neuron
% is defined by set of three cone weights and a 1-D Gabor spatial weighting
% function (time is not explictly simulated). Within each pair of linear
% neurons, cone weights are identical, and RFs are in quadrature. The 
% standard energy model of complex cells uses a pair of linear neurons 
% that are responsive to achromatic modulations. The innovation of the color
% energy model is the inclusion of a second pair of linear neurons, with 
% quadrature-phase RFs, that can be tuned for a different color direction.
%
% The complex cell combines inputs according to an ENERGY calculation
% (a^2+b^2+c^2+d^2) or an NONLINEARENERGY calculation ((a+c)^2+(b+d)^2)).
% This is identical to the energy calculation plus two terms: 2ac+2bd 
% that boost the response if the achromatic and chromatic cosine components
% (a,c) have the same sign or if the achromatic and chromatic sine components
% (b,d) have the same sign.
%
% If the NAKARUSHTON flag is set to true in nested function, SimResp, the 
% result of the energy calculation is passed through a Naka-Rushton function
% before being returned. This feature is included to show that the 
% isoresponse method is immune to static nonlinearities. 
%
% Stimulus:
% The stimulus is a compound grating consisting of two spatial modulations
% (e.g. L+M+S and L-M) of identical spatial frequency, fixed relative
% contrast, and parametrically varied relative phase.
% On any given trial, the contrast of the two components, and their relative 
% phase, is fixed. Across trials, the overall contrast is varied to
% find the contrast that produces the target response. 
%
% Experimental procedure:
% First, each stimulus component is presented in isolation and the contrast
% is titrated to obtain the user-defined target firing rate (CRITFR).
% These relative contrasts are fixed for the remainder of the experiment.
%
% Second, the two components from phase 1 are added in 0 and 180° phase 
% to create two compound gratings (e.g. a bright-R/dark-G grating and a 
% bright-G/dark-R grating).
% The overall contrast of each compound grating is adjusted to find 
% the contrast that produces the target firing rate. This contrast need not
% be identical in both directions. For example, an LN neuron would fail to
% respond to one of these stimuli (stimulus modulations would cancel in 
% one of the two phases). 
%
% Third, N phases are selected and isoresponse measurements are made for
% stimuli with those phases. For an energy model cell, the isoresponse points
% lie on a circle. For a linear neuron the line on lines (either horizontal
% or vertical. For a neuron that is specialized for encoding in-phase
% modulations of luminance and chromatic signals, the isoresponse contour
% will be closer to the origin for the 0 and 180° phases than the 90 and
% 270° phases.
%
% Note, the predictions from the energy model hold for any
% combination of cone weights and stimuli, provided the projections of 
% the stimuli span the plane defined by the two sets of cone weights.
% We do not need to assume that the neurons is tuned for and particular
% color directions.

coneweights = [3 1 0; 10 -10 0 ]; % cone weights of two linear mechanisms that combine to create a color-complex cell
stimcc = [1 1 1;.1 -.1 0]; % two stimuli parametrically shifted in relative phase
%stimcc = unifrnd(-1,1,2,3);
%stimcc*coneweights' % checking how strongly each stimulus acts on each mechanism


if rank(coneweights') == 1
   disp('neuron is only sensitive to one color direction')
end

whichmodel = 'NONLINEARENERGY';
%whichmodel = 'ENERGY';
GAMUTEDGE = 100;
CRITFR = 20;

PIXELS = linspace(0,6*pi,100);
RFtemplate = sin(PIXELS).*normpdf(PIXELS,mean(PIXELS),2); 
RFtemplate(2,:) = cos(PIXELS).*normpdf(PIXELS,mean(PIXELS),2);

RFs = zeros(size(coneweights,2),length(PIXELS), 2, size(coneweights,1)); % cones x space x sin/cos x mechanism
for i = 1:size(coneweights,1) % mechanisms
    for j = 1:size(coneweights,2) % cones
        for k = 1:size(RFtemplate,1) % sin/cos
            RFs(j,:,k,i) = RFtemplate(k,:).*coneweights(i,j);
        end
    end
end

% Simulated experiment
% First, finding the contrast of the two components that produce the
% criterion firing rate.
componentcontrasts = [nan nan];

for whichcomponent = 1:2
    componentcontrasts(whichcomponent) = findContrast([stimcc(whichcomponent,:); 0 0 0], [0 0], PIXELS, RFs, whichmodel, CRITFR, GAMUTEDGE);
end
if any(componentcontrasts == GAMUTEDGE)
    disp('One of the components does not drive a criterion response');
    keyboard
end
stimcc = diag(componentcontrasts)*stimcc; % Destructively modifying stimcc so that both components produce same response

% Now finding contrast of the compound grating in the 0 and pi phases that
% produce the criterion response.
phases = [0 pi];
startingcontrasts = [nan nan];
for phase = phases
    startingcontrasts(phases == phase) = findContrast(stimcc, [0 phase], PIXELS, RFs, whichmodel, CRITFR, GAMUTEDGE);
end

% Figuring out which phases to test based on the relative sensitivity to
% the 0° and 180° phase-shifted stimuli.
phases = linspace(0,pi,21);
[x,y] = pol2cart(phases,1); % get a bunch of points on the unit circle
x = x.*startingcontrasts(1); y = y.*startingcontrasts(2);
phases = atan2(y,x); phases = phases*2; phases(phases<0) = 2*pi+phases(phases<0);

% Now the actual experiment
isorespcontrasts = zeros(size(phases));
for phase = phases
    isorespcontrasts(phases == phase) = findContrast(stimcc, [0 phase], PIXELS, RFs, whichmodel, CRITFR, GAMUTEDGE);
end
figure; subplot(2,1,1); 
LOOG = isorespcontrasts >= GAMUTEDGE;
polar(phases(~LOOG),isorespcontrasts(~LOOG),'k-o'); hold on;
polar(phases(LOOG),isorespcontrasts(LOOG),'ro');

subplot(2,1,2); hold on;
[x,y] = pol2cart(phases/2, isorespcontrasts);
plot(x(~LOOG)./max(x),y(~LOOG)./max(y),'k-o');
plot(x(LOOG)./max(x),y(LOOG)./max(y),'ro');
axis equal


% Nested functions
function c = findContrast(stimcc, phases, PIXELS, RFs, whichmodel, criterion, GAMUTEDGE)
   
    try
         x0 = 0;
         [c,fval,exitflag] = fminsearch(@nestedfun, x0); 
         if (c == x0 & exitflag == 1) | c > GAMUTEDGE
             c = GAMUTEDGE;
         end
     catch
         keyboard
     end
     
     function err = nestedfun(x0)
         err = (SimResp(MkStim(x0, stimcc, phases, PIXELS), RFs, whichmodel)-criterion).^2;
     end
end

function stim = MkStim(contrast, stimcc, phases, x)
     stim = zeros(size(stimcc,2),length(x),size(stimcc,1));
     stim(:,:,1) = contrast*stimcc(1,:)'*cos(x+phases(1));
     stim(:,:,2) = contrast*stimcc(2,:)'*cos(x+phases(2));
     stim = sum(stim,3);
end

function response = SimResp(stim, RFs, whichmodel)
    NAKARUSHTON = true;
    tmp =[];
    for k = 1:size(RFs,3) % sin/cos
        for l = 1:size(RFs,4) % mechanisms
            tmp = [tmp, sum(sum(RFs(:,:,k,l).*stim))];
        end
    end
    if strcmp(whichmodel, 'ENERGY')
        response = sqrt(sum(tmp.^2));
    elseif strcmp(whichmodel, 'NONLINEARENERGY')
        response = (abs(tmp(1))+abs(tmp(2)))^2+(abs(tmp(3))+abs(tmp(4)))^2;
    else
        error('Unknown model type');
    end
    if NAKARUSHTON
        response = 200*response/(response+50);
    end
end



