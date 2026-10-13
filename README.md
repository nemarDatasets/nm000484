[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000484-blue)](https://doi.org/10.82901/nemar.nm000484)

# Example intracranial EEG with qHFO detections (Stacey lab, Mendeley Data bt58z3yhdy)

The intracranial example recording distributed with the code of Ren, Gliske, Brang and Stacey (2019), "Redaction of false
high frequency oscillations due to muscle artifact improves specificity to epileptic tissue" (Clin Neurophysiol): 61 intracranial
channels, 3000001 samples at 5000 Hz (10.0 min), with 357 automated qHFO
detections (start/stop sample and channel). The release uses it to demonstrate the wavelet-entropy feature that separates
EMG-contaminated false HFOs from real HFOs (threshold 8.956 in the release script).

## Source
- Mendeley Data doi:10.17632/bt58z3yhdy.1 (contributor William Stacey), licence CC0 1.0. `examples.zip` = `intracranialEEG.mat`
  (MATLAB v7.3, `data` int16 61 x 3000001), `qHFOdetection.mat` (startIdx, stopIdx, chanIdx), `data.mat` (19-channel scalp EEG
  example for the scalp EMG detector), MATLAB code and `Readme.docx`. All kept byte-identical in `sourcedata/mendeley_bt58z3yhdy/`.
- Sampling frequency: `fs = 5000` in `example_intracranial_feature.m` ("sampling frequency of raw iEEG"); the scalp example
  uses `FS = 5000` as well.

## Content
- `sub-01/ieeg/sub-01_task-rest_ieeg.*`: the 61 intracranial channels as BrainVision IEEE float32, integer values unchanged.
  The release gives no channel names (channels are `ch01`..`ch61` = matrix rows), no electrode types or positions, no reference
  and no unit (`units` = n/a; the amplitudes, std ~250, are consistent with microvolts but this is not stated). Channel type OTHER.
- `events.tsv`: one row per qHFO detection (onset/duration from the release sample indices, channel). Whether a detection is
  an EMG-false HFO is computed by the release code (not released as a table).
- The 19-channel scalp EEG example (`data.mat`) is in sourcedata only (no channel names or montage released).

## Participants
One patient under intracranial monitoring; no demographic or clinical information in the release.

## Privacy
The release holds anonymous matrices and sample indices only; no dates (`acq_time` n/a), names or identifiers.
