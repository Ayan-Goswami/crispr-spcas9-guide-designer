# CRISPR SpCas9 Guide RNA Designer

A lightweight, zero-dependency bioinformatics tool for designing and validating CRISPR SpCas9 guide RNAs (sgRNAs). Available as both a standalone Python script and an interactive client-side web application.

This project automates the first step of CRISPR gene editing: scanning a target DNA sequence for Cas9 PAM sites (`NGG`), extracting the 20-nucleotide upstream spacer, and filtering candidates by thermodynamic stability (GC content).

## Features

* **Double-Stranded Scanning:** Identifies `NGG` PAM sites on the forward strand and `CCN` PAM sites on the reverse complement strand.
* **Thermodynamic Filtering:** Calculates GC content for every 20-nt spacer and filters for the optimal 40%–60% window to prevent weak binding or hairpin formation.
* **Smart Ranking:** Sorts viable guides by proximity to the ideal 50% GC biochemical sweet spot.
* **Zero Dependencies:** Pure Python implementation. No external libraries, databases, or API keys required.
* **Web App Interface:** Includes a single-file HTML/JS web application with a sequence inspector and CSV export.

## Repository Structure

```text
crispr-spcas9-guide-designer/
├── index.html               # Interactive web application (Vanilla HTML/JS/CSS)
├── src/
│   └── crispr_designer.py   # Command-line Python script
└── README.md
