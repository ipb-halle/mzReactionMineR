## Summary
This repository contains an R package for the analysis of mass spectrometry data that was pre-processed using other software (mzmine, MetaboScape, etc.). It helps in integrating and analyzing the processed data efficiently within the R environment. It does so by utilizing SummarizedExperiment objects to manage and manipulate the data.

## Terminology
- rt - retention time as in the time where a feature eluted from the chromatography column.
- mz - the mass-to-charge ratio of a measured feature
- ims - ion mobility spectrometry. An additional separation dimension based on the mobility of ions in a gas phase.
- cc - collision cross-section. A measure of the effective surface area of an ion as it travels through a gas. based on the ims value.
- ms2 - tandem mass spectrometry, a technique where ions are fragmented and the resulting fragments are analyzed to provide structural information about the original ion.
- feature - a detected signal in the mass spectrometry data, typically characterized by its retention time (rt), mass-to-charge ratio (mz), and optionally ion mobility (ims) and collision cross-section (cc).
- spectral similarity - a measure as to how similar two features are. Usually based on their ms2 spectra.
- alignment - the process of matching features across different samples to account for variations in retention time, mass-to-charge ratio, and optionally ion mobility or spectral similarity, ensuring that the same feature is consistently identified across all samples.
- mass difference network - a network where nodes represent features and edges represent the mass differences between them, often used to infer potential biochemical relationships or reactions between the features.
- is - internal standard. A compound added to samples in a known quantity to help with the normalization of the measured features.
- blank - a sample that does not contain the analyte of interest and is used to identify background signals and contaminants.
- mzmine - aka MZmine or mzMine. A popular open-source software for mass spectrometry data processing, including feature detection, alignment, and quantification.
- MetaboScape - a software developed by Bruker for the analysis of mass spectrometry data, including feature detection, alignment, and quantification, often used in conjunction with Bruker instruments.
- SIRIUS - a software for the analysis of mass spectrometry data, particularly for the identification of molecular formulas and structures based on tandem mass spectrometry (MS2) data.
- feature table - a structured representation of detected features in the mass spectrometry data. Rows typically correspond to features. Columns contain various information about the feature and intensities across different samples.

## Architecture
Contains mainly helper functions and wrappers that facilitate the manipulation, analysis, and visualization of mass spectrometry data within the R environment. Any utility functions that should not be used by the end-user directly are written in the utils.R file.

## Task planning and problem-solving
- Before each task, you must first complete the following steps:
  1. Provide a full plan of your changes.
  2. Provide a list of behaviors that you'll change.
  3. Provide a list of test cases to add.
- Before you add any code, always check if you can just re-use or re-configure any existing code to achieve the result.
- Use existing utility functions whenever possible to avoid code duplication and maintain consistency across the project.
- Create utility functions for common tasks to promote code reuse and maintainability.
- Utilize functions from other packages and libraries whenever possible, if they are well maintained and suitable for the task at hand.

## Coding guidelines
- Always write clear and concise code that is easy to understand and maintain.
- Follow consistent naming conventions for variables, functions, and files across the repository.
- Document all functions with meaningful comments and usage examples.
- Avoid code duplication by reusing existing functions and utilities.
- Write unit tests for all new functionality to ensure correctness and reliability.
- Ensure that the code adheres to the overall architecture and design principles of the project.
- add error handling to ensure that the code gracefully handles unexpected inputs and edge cases.

## Tooling

When running Git commands, use:

C:\Program Files\Git\bin\git.exe

When running R scripts, use:

C:\Program Files\R\R-4.6.1\bin\Rscript.exe

Do not assume these executables are available via PATH.
Always use the full paths above when executing commands.