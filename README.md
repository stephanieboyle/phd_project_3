# Why do high notes feel small?

### PhD project 3 · Auditory pitch and visual size

People often pair high-pitched sounds with small objects and low-pitched sounds with large ones. This project asks **when those cross-sensory associations show up in brain activity**, and how they relate to response times.

**[Read the preprint](https://www.biorxiv.org/content/10.1101/423939v1.abstract)** · **[Explore the data on OSF](https://osf.io/e2nv5/)**

## What I investigated

A modified implicit association task paired sounds and circles with congruent or incongruent response mappings. Behavioural measures and 128-channel EEG were analysed using two single-trial approaches, including mutual information analysis.

## What the study found

Participants responded faster with congruent pairings. EEG effects related to pitch and visual size emerged early, and auditory effects predicted trial-by-trial response times. Findings were consistent across the two EEG analysis approaches, supporting an early perceptual contribution to these associations. See the preprint for the full methods and results.

## What's in this repository?

| Files | Purpose |
|---|---|
| [Manual preprocessing](analysis%20scripts/AnalysisP2_Preprocessing_Manual.m) | Interactive EEG preprocessing steps. |
| [Automatic preprocessing](analysis%20scripts/AnalysisP2_Preprocessing_Automatic.m) | Automated preprocessing stages. |
| [Behavioural and EEG analysis](analysis%20scripts/AnalysisP2_CircleExperiment_NEW.m) | Analysis for the pitch–size task. |
| [Mutual information](analysis%20scripts/AnalysisP2_NewMI.m) | Single-trial information analysis using GCMI functions. |
| [Paper figures](analysis%20scripts/AnalysisP2_PaperFigures.m) | Figure-generation analysis. |
| [Original ZIP](analysis%20scripts.zip) | Archived scripts bundle. |

Some original filenames use “P2”; they are retained as part of this project's original source.

## Using the code

This is an archive of the original research code. The scripts contain paths to the original lab drives and expect study data and intermediate MATLAB files outside this repository. Start with the linked data and research write-up, then review paths and dependencies before running individual sections.

The source uses MATLAB and FieldTrip functions alongside lab-specific helpers. The original MATLAB/toolbox versions are not pinned here, and the analysis has not been rerun as part of this README refresh. The original scripts and ZIP archives are preserved.

The mutual-information script also references the GCMI toolbox and additional helper functions.

## Related research

[Project 1](https://github.com/stephanieboyle/phd_project_1) · [Project 2](https://github.com/stephanieboyle/phd_project_2)

[More about me](https://stephanieboyle.github.io/)
