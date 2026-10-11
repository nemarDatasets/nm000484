"""o5c lane D: src-literature-83 (Mendeley bt58z3yhdy, Stacey lab, 'Redaction of false HFOs due to muscle artifact', CC0):
example intracranial EEG (intracranialEEG.mat 'data' [3000001 x 61] int16, fs = 5000 Hz per example_intracranial_feature.m)
with qHFO detections (qHFOdetection.mat startIdx/stopIdx/chanIdx) -> one BIDS-iEEG run (BrainVision float32 of the
int16 values, unit not stated in the release -> n/a) + events.tsv (one event per detection). The scalp example
(data.mat, 19 ch) and the MATLAB code stay in sourcedata/. Usage: python o5c_D_convert_l83.py <record root> <bids out>"""
import json, os, shutil, sys
import numpy as np
import h5py
import scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b3w10_common import wtsv, wjson, write_run, sha256_file

ROOT, OUT = sys.argv[1], sys.argv[2]
D = os.path.join(ROOT, "_x", "examples.zip.d")
SF = 5000.0
sub, stem = "sub-01", "sub-01_task-rest"
od = os.path.join(OUT, sub, "ieeg"); os.makedirs(od, exist_ok=True)
REP = {"sha256_intracranialEEG": sha256_file(os.path.join(D, "intracranialEEG.mat")), "sha256_qHFO": sha256_file(os.path.join(D, "qHFOdetection.mat"))}
with h5py.File(os.path.join(D, "intracranialEEG.mat"), "r") as f:
    d = f["data"]
    REP["h5_shape"] = list(d.shape); REP["h5_dtype"] = str(d.dtype)
    X = np.asarray(d[()])  # stored [3000001 x 61] (samples x channels) in the HDF5 view
if X.shape[0] < X.shape[1]:
    X = X.T
n, nch = X.shape
assert nch == 61, X.shape
names = [f"ch{i+1:02d}" for i in range(nch)]
hfo = sio.loadmat(os.path.join(D, "qHFOdetection.mat"), squeeze_me=True)
st = np.asarray(hfo["startIdx"]).astype(np.int64).ravel(); sp = np.asarray(hfo["stopIdx"]).astype(np.int64).ravel(); ch = np.asarray(hfo["chanIdx"]).astype(np.int64).ravel()
markers = [("Stimulus", f"qHFO_ch{int(c):02d}", int(s - 1), max(1, int(e - s + 1))) for s, e, c in zip(st, sp, ch)]
ch_extra = [{"low_cutoff": "n/a", "high_cutoff": "n/a", "reference": "n/a", "group": "n/a", "status": "good",
             "description": "intracranial channel %d of the release matrix (row index; the release gives no channel names, locations or types)" % (i + 1)} for i in range(nch)]
sidecar = {"TaskName": "rest", "TaskDescription": "No task is described in the release (example intracranial EEG recording used to demonstrate the wavelet-entropy feature that separates EMG-false HFOs from real HFOs)",
           "Manufacturer": "n/a", "PowerLineFrequency": 60, "SoftwareFilters": "n/a", "HardwareFilters": "n/a", "iEEGReference": "n/a", "ElectrodeManufacturer": "n/a",
           "ElectricalStimulation": False, "RecordingType": "continuous", "InstitutionName": "University of Michigan (Stacey lab; release contributor William Stacey)",
           "Notes": "Release intracranialEEG.mat 'data' int16 [3000001 x 61]; sampling frequency 5000 Hz from the release script example_intracranial_feature.m ('fs = 5000; %% sampling frequency of raw iEEG'). The unit of the integer values is not stated (units = n/a); values written unchanged as float32. The release does not say whether the 61 channels are ECoG or SEEG contacts, so the channel type is OTHER (intracranial, unspecified)."}
r = write_run(od, stem, X.astype(np.float32), SF, names, ["n/a"] * nch, ["OTHER"] * nch, sidecar, ch_extra=ch_extra, markers=markers,
              comment="o5c_D_convert_l83.py: intracranialEEG.mat data int16 -> IEEE float32 (values unchanged)")
REP["run"] = r
# events: qHFO detections
wtsv(os.path.join(od, f"{stem}_events.tsv"), ["onset", "duration", "trial_type", "channel", "sample", "stop_sample"],
     [[round((s - 1) / SF, 6), round((e - s + 1) / SF, 6), "qHFO", names[int(c) - 1], int(s - 1), int(e - 1)] for s, e, c in zip(st, sp, ch)])
wjson(os.path.join(OUT, "task-rest_events.json"), {"trial_type": {"Description": "qHFO detection from the release qHFOdetection.mat (automated high-frequency-oscillation detector output; some are EMG-false HFOs, see the release code and the article)", "Levels": {"qHFO": "detected HFO candidate"}},
                                                   "channel": {"Description": "Channel of the detection (release chanIdx, 1-based row -> chNN)"}, "sample": {"Description": "Release startIdx (0-based here)"}, "stop_sample": {"Description": "Release stopIdx (0-based here)"}})
REP["n_events"] = int(len(st))
wtsv(os.path.join(od, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group"], [[nm, "n/a", "n/a", "n/a", "n/a", "n/a"] for nm in names])
wjson(os.path.join(od, f"{sub}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a", "iEEGCoordinateSystemDescription": "No electrode names, positions or types are released; one entry per matrix row with x, y, z = n/a."})
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "group"], [[sub, "n/a", "n/a", "epilepsy (intracranial monitoring; release gives no details)"]])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "Not given in the release", "Units": "year"}, "sex": {"Description": "Not given in the release"}, "group": {"Description": "Patient with intracranial EEG monitoring (example recording of the release)"}})
SD = os.path.join(OUT, "sourcedata", "mendeley_bt58z3yhdy"); os.makedirs(SD, exist_ok=True)
for f in os.listdir(D):
    p = os.path.join(D, f)
    if os.path.isfile(p):
        shutil.copyfile(p, os.path.join(SD, f))
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
wjson(os.path.join(OUT, "code", "conversion_report.json"), REP)
shutil.copyfile(os.path.abspath(__file__), os.path.join(OUT, "code", os.path.basename(__file__)))
print(json.dumps(REP, default=str)[:1500]); print("CONVERT_DONE")
