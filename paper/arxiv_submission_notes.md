# arXiv Submission Notes for "The Arrow of Time from Petz Recovery"

## Metadata

**Title:** The Arrow of Time from Petz Recovery

**Authors:** Sheng-Kai Huang

**Email:** akai@fawstudio.com

**Affiliation:** Independent Researcher

## arXiv Categories

- **Primary:** quant-ph (Quantum Physics)
- **Cross-list:** math-ph (Mathematical Physics), cond-mat.stat-mech (Statistical Mechanics)

## Abstract

The Petz recovery map is simultaneously the unique Bayesian retrodiction
functor, a near-optimal quantum error correction decoder, and the reverse
channel in quantum Crooks theorems. We show that these roles are not
coincidental: the conditions for exact Petz recovery, vanishing entropy
production, the quantum Markov condition, and the quantum eraser are all
equivalent, forming a four-way equivalence chain (Eraser <-> Retrodiction
<-> QEC <-> Thermodynamics) unified by the Petz map. Writing the recovery
infidelity as tau = 1 - F, we prove two new results: (1) the bound
F^2 >= exp(-Delta D) of Junge et al. is exactly saturated if and only if
the channel acts as a quantum sufficient statistic (Theorem 1), and
(2) sqrt(tau) satisfies a triangle inequality under channel
composition---a structural property unique to the standard Petz map that
fails for the rotated variant (Theorem 2). The equivalence chain, together
with these theorems, provides a dictionary that lets results in any one of
the four domains be translated into the other three.

## Comments Field

6 pages, 1 table; Supplemental Material included (20 sections)

## Files to Submit

1. `petz_recovery_unification.tex` (main paper)
2. `petz_recovery_supplemental.tex` (Supplemental Material)
3. `petz_recovery_supplemental.pdf` (compiled SM, if submitting as ancillary)

Note: The main paper uses inline \thebibliography (no .bib file needed).
The Supplemental Material also uses inline \thebibliography.
Compile with: pdflatex (three passes, no bibtex needed for main paper).

## MSC Codes

- 81P47 (Quantum information, communication, networks)
- 94A40 (Channel models in information theory)
- 81P45 (Quantum information, open systems, decoherence)

## PACS Numbers (if requested)

- 03.67.Pp (Quantum error correction)
- 03.65.Ta (Foundations of quantum mechanics; measurement theory)
- 05.70.Ln (Nonequilibrium and irreversible thermodynamics)

## Journal Target

Physical Review Letters (PRL) -- revtex4-2 format, 6 pages

## Key References to Highlight

- Petz (1986, 1988): Original recovery map
- Parzygnat-Buscemi (2023): Retrodiction uniqueness
- Junge-Renner-Sutter-Wilde-Winter (2018): Universal recovery bound
- Fawzi-Renner (2015): Quantum conditional mutual information
- Kwon-Kim (2019): Quantum Crooks via Petz
- Bai-Buscemi-Scarani (2024): Concurrent work on stochastic entropy production

## Pre-submission Checklist

- [x] Paper compiles without errors (6 pages)
- [x] No undefined citations
- [x] No overfull hboxes
- [x] All references verified (journals, volumes, pages, years)
- [x] Supplemental Material sections S11, S17, S18, S19, S20 all present
- [x] Concurrent work (Buscemi et al. 2024) acknowledged in Introduction and Sec. III
- [x] Author info complete (name, email, affiliation)
- [x] Crooks (1999) and Acharya (2025) now cited in text
- [x] Note1/Note2 footnote bibitems included for revtex4-2 compatibility
