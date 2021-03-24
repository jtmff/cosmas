% This script shows how to analyse calcium transient (CaT) alternans in a rapidly-paced heart. It
% visualizes and quantifies the discordant alternans present.
clear parameters params

params.extension = '.png';
params.binningFactor = 2;
[imStackCaAlt, atCaAlt] = COSMAS.readFolder('data/alternansDiscordantCa_bcl85', params);

parameters.spikesPointDown = false;
parameters.baselineDefinition = 'first'; % if we want to measure alternans in amplitude of CaT upstroke, we measure baseline as the first element in each segmented calcium transient.
[baselineCaAlt, amplitudeCaAlt, durationCaAlt, activationMapsCaAlt] = COSMAS.analyseRegularPacing(imStackCaAlt, 85, parameters);


alternansMap = COSMAS.getAlternans(amplitudeCaAlt.maps, 'sMAPE');

% We plot the alternans map - cold colours correspond to nodal lines, where
% no alternans is present
figure(1); clf
imagesc(alternansMap); colorbar
title(['alternans map, mean alternans = ' num2str(mean(alternansMap(:)))]);

% As an additional visualization, we show the mean CaT for even/odd wave
% passes - where these differ, alternans is present
figure(2); 
title('CaT amplitude for even/odd waves');
subplot(1,2,1);
imagesc(amplitudeCaAlt.mapMeanEven); colorbar; caxis([30, 170]);
subplot(1,2,2);
imagesc(amplitudeCaAlt.mapMeanOdd); colorbar; caxis([30, 170]);

% Note: average alternans should be ideally taken as an average of an
% alternans map. An alternative and potentially appealing way is:
alternansData = COSMAS.getAlternans(amplitudeCaAlt.data, 'sMAPE');
%, where alternans on average CaT per-wave-pass is measured. However, this
% is confounded by the discordant nature of alternans; just
% taking averages between even and odd averages for calcium transient is
% overshadowing the discordance. If unclear, imagine a discordant alternans
% where the left half of recording is in the opposite phase of the right
% half (and the sides are otherwise identical) - then there is no
% difference betwen beats in the calcium transient amplitude when averaged
% over the field of view, yet there can be pronounced alternans).
% - However, one may in fact leverage this "problem" for high-throughput
% detection of discordant (versus concordant) alternans - whenever the
% variable 'alternansData' above differs markedly from the value of
% 'alternansMap' variable, it indicates that alternans is more likely to be
% discordant