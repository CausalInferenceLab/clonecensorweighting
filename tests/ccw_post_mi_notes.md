# CCW simulated example: post-MI antiplatelet initiation window

This note accompanies `ccw_post_mi_example.R`.

## Clinical story

- **Time zero**: discharge after myocardial infarction (MI; heart attack).
- **Exposure strategy 1**: start an antiplatelet medication within 30 days.
- **Exposure strategy 2**: do not start within 30 days.
- **Outcome**: time to a serious cardiovascular event (illustrated as death or urgent cardiovascular readmission) within 180 days.

## Variables in the simulation

- `age`: age in years.
- `bleed`: high bleeding risk / contraindication indicator.
- `frailty`: overall frailty / comorbidity burden.
- `large_mi`: whether the index heart attack was relatively large/severe.
- `L`: daily latent ischemic instability. Think of this as a summary of symptoms, physician concern, ECG/troponin instability, and overall “how worried are we today?”
- `start_day`: the day the medication starts **at the end of the day**.
- `event_day`: the day the outcome occurs **during the day**.

## Why this creates confounding

Patients with higher `L` are:
1. more likely to be started on treatment, and
2. more likely to have the outcome.

Therefore, simple treated-vs-untreated comparisons are biased.

## Why CCW is needed

If we compare “started within 30 days” vs “did not start within 30 days” directly from baseline, the treated group must survive long enough to start treatment. That produces immortal-time bias.

CCW fixes this by:
1. cloning each person into both strategies at day 0,
2. censoring each clone when observed care deviates from the assigned strategy,
3. weighting to account for informative artificial censoring,
4. fitting a weighted outcome model.

## Weighting choices in the script

- `w_ipcw`: standard stabilized inverse probability of censoring weights.
- `w_ipcw_trunc`: truncated version of `w_ipcw` for stability.
- `w_ow`: overlap-style sensitivity analysis based on cumulative uncensoring probability.

`w_ipcw_trunc` is the default recommendation for the main analysis.
