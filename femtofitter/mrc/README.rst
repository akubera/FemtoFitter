==============================
Momentum Resolution Correction
==============================


There are two momentum resolution correction procedures:


``mrc::CorrFctnRatio``

    This is classic momentum resolution correction solution.

    Objects are created with 2 pairs of histograms: the correlation function
    created using the montecarlo "generated" momentum, and the "reconstructed"
    momentum.

    By making the ratio of the generated correlation function to the reconstructed
    correlation function, we create the "unsmearing factor", which is multiplied
    to the data to correct.

    Note that this requires the use of a femto-weight generator to ensure that
    femto effects are taken into account by the two correlation functions, leaving
    only detector effects in the correction factor.


``mrc::MomentumMap``

    The momentum-map method does not create correlation functions, but rather
    creates a map of the reconstructed momentum space to the generated
    momentum space.
    For 1D Correlation functions, this is simply a 2D histogram of $qinv_{rec}$
    vs :math:`qinv_{gen}`. The 3D case uses a 6D histogram to map each bin of
    reconstructed q_{o,s,l} to the generated q_{o,s,l}.
    This is efficiently handled with the THnSparse histogram class.

    An additional step of normalization is required before using the momentum-map
    correction.
    We assume that the reconstructed values come from a distribution of the
    generated "source" momentum.
    With infinite detector efficiency, we would expect this distribution to be
    a delta function mapping each bin in the generated space to the same bin
    in reconstructed space.
    In reality, the distribution is near-exponential, centered at the corresponding
    bin in generated-space.
    We simply normalize this distribution to 1, allowing the pair-count in each
    data bin to be multiplied by this distribution, correcting the detector
    effects.
    The final result is the sum of the smearing (or rather, unsmearing) of each
    of the reconstructed (smeared) bins into the generated (unsmeared) space.
    The integral of the corrected histogram should have the same number of
    entries as the original, though this method *will* push some of entries at the
    edges into the overflow bins.

    The corrected histograms are produced by multiplying the counts in each bin
    of the correlation function numerator and denominator histograms by the
    corresponding
    the valuesapt list --upgradableapt list --upgradable in
    the values of this this momentum-map This should
