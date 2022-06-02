% AMT - Publication-specific models
%
%  Peripheral models
%     DAU1996               - Linear filtering for monaural masking (basic)
%     DAU1997               - Linear filtering for monaural masking (improved)
%     HOHMANN2002           - Invertible Gammatone filterbank
%     LOPEZPOVEDA2001       - Dual resonance non-linear (DRNL) filterbank
%     LYON2011              - Cascade of asymmetric resonators with fast-acting compression (CARFAC) model
%     VERHULST2012          - Cochlear transmission-line model (basic)
%     VERHULST2015          - Cochlear transmission-line model (improved)
%     VERHULST2018          - Cochlear transmission-line model (improved, incl. brainstem)
%     ZILANY2007            - Auditory-nerve filterbank (basic)
%     ZILANY2014            - Auditory-nerve filterbank (improved)
%     BRUCE2018             - Auditory-nerve filterbank (improved synapse) 
%
%   Temporal-modulation sensitivity
%     CARNEY2015            - Brainstem processing
%     ROENNE2012            - Simulate auditory brainstem responses (ABRs)
%     EWERT2000             - Modulation filterbank (based on EPSM)
%     KING2019              - Modulation filterbank (based on nonlinear processing)
%     RELANOIBORRA2019      - Modulation filterbank (based on DRNL)
%
%   Binaural processing
%     LINDEMANN1986         - Binaural activity map based on cross-correlation
%     BREEBAART2001         - Binaural masking level differences
%     TAKANEN2013           - Binaural count-comparison model
%
%   Monaural speech perception
%     JOERGENSEN2011         - Speech-based envelope power spectrum (EPSM)
%     JOERGENSEN2013         - Speech-based envelope power spectrum (multi-resolution EPSM)
%     TAAL2011               - Short-time objective intelligibility
%
%   Binaural speech perception
%     CULLING2004            - Binaural speech intelligibility
%     HAUTH2020              - Blind equalization cancellation model
%     JELFS2011              - Binaural speech advantage
%     LECLERE2015            - Binaural useful-to-detrimental ratio for a reverberated target
%     LAVANDIER2022          - Compute the binaural 'effective' target-to-interferer ratio
%     PRUDHOMME2020          - Compute the effective SNR taking into account harmonic cancellation
%     VICENTE2020            - Computes the effective SNR taking into account audibility (audiogram)
%     VICENTE2020NH          - Compute the effective SNR taking into account BU and BE
%
%   Perceptual similarity
%     OSSES2021              - Monaural perceptual similarity
%     MCKENZIE2021           - Binaural perceptual similarity
%
%   Loudness
%     MOORE1997              - Loudness model for stationary signals
%     GLASBERG2002           - Loudness model for time-variant signals
%     CHEN2011               - Fast excitation pattern estimation
%     MOORE2016              - Binaural loudness model
%
%   Spatial perception
%     DIETZ2011             - Sound lateral direction
%     MAY2011               - Azimuthal localization of concurrent talkers
%     KELVASA2015           - Azimuthal localization in cochlear-implant listeners
%     LANGENDIJK2002        - Median-plane localization probability
%     BAUMGARTNER2013       - Localization in saggital planes (simple)
%     BAUMGARTNER2014       - Localization in saggital planes (robust, linear periphery)
%     BAUMGARTNER2016       - Localization in sagittal planes (robust, nonlinear periphery)
%     HASSAGER2016          - Sound externalization based on interaural level differences
%     BAUMGARTNER2017       - Sound externalization model based on monaural spectral cues
%     BAUMGARTNER2021       - Sound externalization based on perceptual-based decisions
%     LI2020                - Sound externalization in reverberant spaces
%     LLADO2022             - Neural network localization
%     GEORGANTI2013         - Distance estimation
%     REIJNIERS2014         - Bayesian spherical sound localization model (basic)
%     BARUMERLI2021         - Bayesian spherical sound localization model (multi-feature)
%     MCLACHLAN2021         - Bayesian dynamic sound localization model
%     WIERSTORF2013         - Sound localization in wavefield synthesis
%     ZIEGELWANGER2013      - Direction-continuous model of time-of-arrival (TOA) in HRTFs (simple)
%     ZIEGELWANGER2014      - Direction-continuous model of time-of-arrival (TOA) in HRTFs (robust)
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/models/Contents.php

% Copyright (C) 2009-2021 Piotr Majdak, Clara Hollomey, and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 1.1.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


