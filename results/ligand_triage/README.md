# Ligand Triage Summary

Descriptor and structural-alert triage for the MD candidate ligands.
This table is not a validated endpoint-specific ADMET prediction; CYP inhibition, hERG, DILI, clearance, and transporter liability still need dedicated prediction or experimental follow-up.

| target | compound_id | compound_name | mw_Da | clogp | tpsa_A2 | hbd | hba | rotatable_bonds | qed | esol_class | lipinski_pass | veber_pass | pains_alert_count | brenk_alert_count | rough_high_gi_absorption | rough_bbb_permeable | triage_label |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CYP1B1 | Mol_11315 | pseudarflavone A | 322.27 | 2.98 | 100.88 | 2 | 6 | 1 | 0.413 | poorly soluble | yes | yes | 0 | 2 | yes | no | acceptable: review alerts/properties |
| SGLT2 | Mol_13144 | justisecundoside B | 530.48 | 0.12 | 162.6 | 4 | 12 | 6 | 0.295 | moderately soluble | no | no | 0 | 1 | no | no | caution: oral rule violation |
| MAO-B | Mol_14056 | chaetocochin C | 305.16 | 4.65 | 26.3 | 0 | 2 | 1 | 0.717 | poorly soluble | yes | yes | 0 | 1 | yes | yes | acceptable: review alerts/properties |

Files:

- `ligand_triage_summary.csv`: compact table for manuscript/supporting summaries.
- `ligand_triage_full.csv`: full descriptor/rule/alert table including SMILES and alert names.
