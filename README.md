# Breast cancer genomic landscape
### Code repository of publication: *Complex rearrangements fuel ER+ and HER2+ breast tumours (Houlahan, Mangiante, Sotomayor-Vivas, Adimoelja et al., 2025)*

## 📂 Repository content ###
This GitHub repository contains the following folders:
```
└─ analysis/: directory containing code used for data analysis
└─ data/: directory containing input data and source data to reproduce figures
|  └─ input_data/: input data
|  └─ source_data/: source data to reproduce figures
|     └─ extended_data/: source data for extended data figs
|     └─ figures/: source data for main figs
|     └─ supplementary_figure/: source data for supplementary figs
└─ extended_data/: directory code to reproduce extended data figs
└─ figures/: directory code to reproduce main figs
└─ supplementary_tables/: directory containing code to reproduce supplementary tables
└─ supplementary_figures/: directory containing R and other code to reproduce supplementary figures
```

## 💻 ENiClust (Ensemble Integrative Clustering)
The source code for IC-subtypes prediction from whole-exome sequencing (WES) or genome sequencing (WGS) is available through [GitHub repo](https://github.com/cancersysbio/ENiClust) and as a [Docker image](https://hub.docker.com/repository/docker/cancersysbio/eniclust-public/general).

Run ENiClust using Docker:
```bash
docker run cancersysbio/eniclust-public:v2.1.0 ENiClust.py --help
```

Run ENiClust using Singularity:
```bash
singularity build eniclust.sif docker://cancersysbio/eniclust-public:v2.1.0

singularity exec eniclust.sif ENiClust.py --help
```

## ✔️ Citation
Kathleen E. Houlahan*, Lise Mangiante*, Cristina Sotomayor-Vivas*, Alvina Adimoelja*, Seongyeol Park, Aziz Khan, Sophia J. Pribus, Zhicheng Ma, Jennifer Caswell-Jin & Christina Curtis. **Complex rearrangements fuel ER+ and HER2+ breast tumours**, *Nature* (2025); doi: [https://doi.org/10.1038/s41586-024-08377-x](https://doi.org/10.1038/s41586-024-08377-x)

## 📨 Contact 
For any questions, please contact Kathleen Houlahan <khoulaha@stanford.edu> or Lise Mangiante <lisem@stanford.edu> or Cristina Sotomayor-Vivas <csotomv@stanford.edu> or Alvina Adimoelja <alvinaa@stanford.edu>
