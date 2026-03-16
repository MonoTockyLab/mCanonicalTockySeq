# mCanonicalTockySeq: Multidimensional Canonical Tocky Analysis of Temporal and Developmental Structure in Single-Cell Transcriptomes

[![Documentation](https://img.shields.io/badge/docs-GitHub_Pages-blue.svg)](https://monotockylab.github.io/mCanonicalTockySeq/index.html)

[![DOI](https://zenodo.org/badge/1183084801.svg)](https://doi.org/10.5281/zenodo.19054871)

<a href="https://monotockylab.github.io/mCanonicalTockySeq/index.html">
<img src="man/figures/mCanonicalTockySeqLogo.jpg" align="center" width="100%">
</a>

**Author:** Dr Masahiro Ono  
**Date:** 16 March 2026

---

## Overview

**mCanonicalTockySeq** is an R package for analysing single-cell transcriptomes in systems where **temporal progression** and **developmental progression** occur simultaneously.

The package was developed around the **Nr4a3-Tocky** system, in which a Fluorescent Timer reporter provides an experimentally anchored readout of strong T-cell receptor (TCR) signalling history.

More generally, the framework constructs a **shared canonical space** constrained by biologically defined landmark populations, and then extracts temporal and developmental structure using experimentally anchored biological readouts.

Rather than assuming that time and development are orthogonal latent variables, **mCanonicalTockySeq** represents them jointly in a biologically informed canonical geometry.

---

## Interactive visualisations

Example interactive visualisations of Nr4a3-Tocky temporal development in the thymus are available at:

**GitHub Pages:**  
[https://MonoTockyLab.github.io/Nr4a3-Tocky-3D-Models/](https://MonoTockyLab.github.io/Nr4a3-Tocky-3D-Models/)


[![Documentation](https://img.shields.io/badge/docs-GitHub_Pages-blue.svg)](https://MonoTockyLab.github.io/Nr4a3-Tocky-3D-Models/)

In the Nr4a3-Tocky thymocyte setting, **mCanonicalTockySeq** enables:

- reconstruction of **Tocky-defined temporal progression**
- simultaneous quantification of **developmental progression toward CD4 and CD8 states**
- identification of **gene-associated developmental trajectories** across Tocky Time
- projection of new datasets, including **cross-species ortholog-translated data**, into an experimentally anchored reference space

---

## Installation

You can install the development version from GitHub:

```r
remotes::install_github("MonoTockyLab/mCanonicalTockySeq")
```

---

## Availability and documentation

* **GitHub repository:** [mCanonicalTockySeq](https://github.com/MonoTockyLab/mCanonicalTockySeq)
* **Package website:** [https://monotockylab.github.io/mCanonicalTockySeq/index.html](https://monotockylab.github.io/mCanonicalTockySeq/index.html)

---

## Citation

If you use `mCanonicalTockySeq` in your research, please cite the associated manuscript and the R package.

### bioRxiv preprint

For a detailed description of the biological applications and underlying methodology of **mCanonicalTockySeq**, please refer to the associated preprint:

> **Canonical Analysis of Fluorescent Timer-Anchored Transcriptomes Resolves Joint Temporal and Developmental Progression**
> Nobuko Irie, Omnia Reda, Yorifumi Satou, Masahiro Ono
> *bioRxiv* (2026). DOI: to be updated

```bibtex
@article{Ono2026mCanonicalTocky,
  author = {Irie, Nobuko and Reda, Omnia and Satou, Yorifumi and Ono, Masahiro},
  title = {Canonical Analysis of Fluorescent Timer-Anchored Transcriptomes Resolves Joint Temporal and Developmental Progression},
  year = {2026},
  journal = {bioRxiv},
  publisher = {Cold Spring Harbor Laboratory},
  doi = {to be updated},
  url = {https://www.biorxiv.org}
}
```

### R package

[![DOI](https://zenodo.org/badge/1183084801.svg)](https://doi.org/10.5281/zenodo.19054871)

```bibtex
@Manual{Ono2026mCanonicalTockySeq,
  title = {mCanonicalTockySeq: Multidimensional Canonical Tocky Analysis of Temporal and Developmental Structure in Single-Cell Transcriptomes},
  author = {Masahiro Ono},
  year = {2026},
  note = {R package version 0.1.0},
  url = {https://github.com/MonoTockyLab/mCanonicalTockySeq}
}
```

<a href="https://monotockylab.github.io/mCanonicalTockySeq/index.html">
<img src="man/figures/mCanonicalTockySeqTransparent.jpg" align="center" width="100%">
</a>

---

## The Ono Lab (MonoTockyLab)

<img src="man/figures/MonoLab.jpg" alt="MonoTockyLab" align="center" width="40%">

**The Masahiro Ono Lab (MonoTockyLab)** develops experimental and computational approaches to study immune cell dynamics, with a particular focus on the temporal regulation of gene expression in T cells.

The lab is known for the development of **Tocky** (*Timer of Cell Kinetics and Activity*), a platform that uses Fluorescent Timer proteins to analyse transcriptional and signalling dynamics *in vivo* at single-cell resolution. Our research integrates mouse genetics, immunology, flow cytometry, single-cell omics, and computational modelling.

Current research directions include:

* cancer immunology and immunotherapy
* temporal mechanisms of T-cell activation, differentiation, and tolerance
* **Foxp3 transcriptional dynamics** and their regulation *in vivo*
* computational methods for time-resolved single-cell analysis, including **CanonicalTockySeq**

**Principal Investigator:** Dr Masahiro Ono

Dr Ono is the creator of **Tocky**, spanning both its transgenic reporter systems and associated analytical frameworks.

---

## Contact and more

**Email**: <a href="mailto:m.ono@imperial.ac.uk"> <img src="https://upload.wikimedia.org/wikipedia/commons/e/ec/Circle-icons-mail.svg" alt="Email" width="10%"> </a>

**Lab homepage**: <a href="http://monotockylab.github.io"> <img src="man/figures/MonoLab.jpg" alt="MonoTockyLab Homepage" align="center" width="30%"/> </a>

**GitHub**: <a href="https://github.com/MonoTockyLab"> <img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png" alt="GitHub" align="center" width="70" height="70"/> </a>

**Twitter / X**: <a href="https://twitter.com/MonoTockyLab"> <img src="https://upload.wikimedia.org/wikipedia/commons/6/6f/Logo_of_Twitter.svg" alt="Twitter" align="center" width="50" height="50"/> </a>

