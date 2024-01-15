# oppni-code-nwc

This repository provides matlab/octave compatible code for running optimized preprocessing and
analytic pipelines in MRI data, including fMRI, ASL, DTI and anatomical data.

- folders are organized around different datatypes, listed as follows:

    oppni-b: For analysis of BOLD fMRI data. Handles both task-based and resting-state at present, with multi-stage QC. 
             Can run and test multiple pipelines with high efficiency. By far the most advanced of the code sets.
    oppni-a: For analysis of arterial spin labelling data. Shares a lot of parallel structure with oppni-b, but includes
             kinectic modelling options for CBF estimation. Currently works for 2D PASL, other sequences are being actively developed.
    oppni-d: For analysis of diffusion-weighted data. Handles single and multi-shell data, with some upwarping capabilities -- may not yet cover all use cases
             Can run a single pipeline, with some hard code-able options. In active development
    oppni-t: For texture-based analysis of T1-weighted structural data etc.
    oppni-x: For containing shared code and group-level analysis scripts

- check each folder for a separate "_readme.txt" file that summarizes current progress & development, as well as versioning for updates.
