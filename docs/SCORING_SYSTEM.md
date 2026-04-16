# AMP Scoring System

AMP-MineFlow uses an 8-point multi-criteria scoring system to identify antimicrobial peptide candidates. Each criterion contributes 1 point (max score = 8).

## Criteria

### 1. Length (1 point)
Peptide length falls within 10–200 amino acids, the typical range for known AMPs.

### 2. Net Charge (1 point)
Net charge at pH 7.0 ≥ +2. Most AMPs are cationic, enabling electrostatic interaction with negatively charged bacterial membranes.

### 3. Hydrophobic Ratio (1 point)
30–70% hydrophobic residues (A, I, L, M, F, V, W). This range balances membrane insertion with aqueous solubility.

### 4. Amphipathic Moment (1 point)
Hydrophobic moment (μH) > 0.25, calculated using the Eisenberg consensus scale with a helical wheel projection (100° rotation). High μH indicates segregation of hydrophobic and hydrophilic faces.

### 5. Disulfide Potential (1 point)
Even number of cysteines (≥ 2), suggesting potential for disulfide bridge formation that stabilizes the peptide structure.

### 6. Sequence Complexity (1 point)
No single amino acid exceeds 30% of total composition, filtering out low-complexity or repetitive sequences.

### 7. Signal Peptide-like N-terminus (1 point)
The first 15 residues contain > 40% hydrophobic amino acids, resembling a signal peptide for secretion.

### 8. Known AMP Motif (1 point)
Sequence contains a motif matching known AMP families (surfactin, iturin, fengycin, subtilin, subtilosin, or plantazolicin).

## Interpretation

| Score | Classification | Recommendation |
|-------|---------------|----------------|
| 7–8   | Strong candidate | Priority for experimental validation |
| 5–6   | Moderate candidate | Consider for secondary screening |
| 4     | Weak candidate | Include in broad screens only |
| < 4   | Not classified | Filtered out |

## References

- Eisenberg D et al. (1982) Nature 299:371-374
- Boman HG (2003) J Intern Med 254:197-215
- Wimley WC & White SH (1996) Nat Struct Biol 3:842-848
