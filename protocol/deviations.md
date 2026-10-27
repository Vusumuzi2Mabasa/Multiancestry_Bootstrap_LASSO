# Prevalidation correction log

During inspection of the initial tuning freeze, the literal scenario name `null`
was found to be interpreted by pandas as its default missing-value token. The
raw per-replicate results and selections were unaffected, but the tuning null-rate
calculation was NaN and incorrectly invoked the documented fallback. No fresh
validation evaluation had run. Validation calibration was stopped, and the initial
tuning and any partial calibration outputs were archived locally under ignored
work/parser_correction_archive. The table reader now retains literal `null` and
asserts that every tuning null rate and objective is finite. Tuning and validation
calibration are rerun under the corrected implementation hash before final
evaluation. The lambda grid, tuning criterion, seed streams, resampling model,
simulation counts and all evaluation metrics are unchanged. The invalid initial
freeze in Git history is superseded by the subsequent corrected freeze.
