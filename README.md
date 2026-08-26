# Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in Human FFPE Tissue

## TABLE OF CONTENTS

- [GENERAL INFORMATION](#general-information)
- [ASSOCIATED PUBLICATION](#associated-publication)
- [RECOMMENDED CITATION](#recommended-citation)
- [USEFUL LINKS](#useful-links)
- [SUPERSPLAT (PREVIEW DATA ON BROWSER)](#3d-gaussian-splats-preview-data-in-browser)
- [ACCESS THE DATASET](#access-the-dataset)
- [FILE ORGANIZATION](#file-organization)
- [REPOSITORY LINKS](#repository-links)
- [FILE LIST](#file-list)
- [ADDITIONAL NOTES/COMMENTS](#additional-notescomments)

## GENERAL INFORMATION

### Title

**Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in Human FFPE Tissue**

### Authors

Alex Y. H. Wong<sup>1,2*</sup>, Yi Daniel Lu<sup>1,2*</sup>, Ziyuan Zhao<sup>1,2,3</sup>, Seo Woo Choi<sup>1,4</sup>, Hojeong Park<sup>5</sup>, Felix Zhou<sup>6</sup>, Zoltan Maliga<sup>1,2</sup>, Yvonne N. A. Anang<sup>1,7</sup>, Soheil R. Talemi<sup>1,2</sup>, Shannon Coy<sup>1,7</sup>, Gaudenz Danuser<sup>6,9</sup>, Sandro Santagata<sup>1,2,3,7</sup>, Clarence Yapp<sup>1,2†</sup> and Peter K. Sorger<sup>1,2,3†</sup>

### Affiliations

<sup>1</sup> Laboratory of Systems Pharmacology, Harvard Medical School, Boston, MA, USA  
<sup>2</sup> Ludwig Centre at Harvard, Harvard Medical School, Boston, MA, USA  
<sup>3</sup> Department of Systems Biology, Harvard Medical School, Boston, MA, USA
<sup>4</sup> LifeCanvas Technologies, Cambridge, Boston, MA, USA 
<sup>5</sup> Broad Institute of MIT and Harvard, Cambridge, MA, USA  
<sup>6</sup> Lyda Hill Department of Bioinformatics, UT Southwestern Medical Center, Dallas, TX, USA  
<sup>7</sup> Department of Pathology, Brigham and Women's Hospital, Boston, MA, USA  
<sup>8</sup> Current address: Tissue Biomarker Laboratory of the Center for Immuno-Oncology, Department of Medical Oncology, Dana-Farber Cancer Institute, Boston, MA, USA  
<sup>9</sup> Current address: Institute for Human Biology, Roche Pharma Research and Early Development (pRED), Roche Innovation Center, Basel, Switzerland



### Author Notes

`*` Authors contributed equally  
`†` Corresponding authors: Clarence Yapp and Peter K. Sorger  
Contact: [peter_sorger@hms.harvard.edu](mailto:peter_sorger@hms.harvard.edu)

### Overview

This repository accompanies a cyclic light-sheet microscopy study for 3D spatial profiling of immunological structures in human specimens. It is intended to organize analysis resources, dataset access details, and links associated with the manuscript.

## ASSOCIATED PUBLICATION

### Publication Details

**Manuscript title:** *Volumetric Cyclic Immunofluorescence for 3D Spatial Profiling of Immune Structures in Human FFPE Tissue*  
**Publication status:** TBD  
**Journal:** TBD  
**Preprint:** [https://doi.org/10.64898/2026.05.17.725158](https://doi.org/10.64898/2026.05.17.725158)   
**DOI:** TBD

## RECOMMENDED CITATION

### Citation Text

Please cite this dataset and repository as:

Wong, A. Y. H. & Lu, Y. et al. (2026). *Cyclic Light-sheet Microscopy for 3D Spatial Profiling of Immunological Structures in Human Specimens*. `{journal / bioRxiv / DOI to be added}`.

## USEFUL LINKS

### External Links

- Publication DOI: `TBD`
- Archived record of this repository: [https://doi.org/10.5281/zenodo.20170967](https://doi.org/10.5281/zenodo.20170967)
- Online data exploration page: `TBD`
- License / restrictions placed on the data: `TBD`

## 3D Gaussian Splats (PREVIEW DATA IN BROWSER)

### Browser Preview

- 20X normal colon: [https://superspl.at/scene/b181b113](https://superspl.at/scene/b181b113)
- virtual H&E: [https://superspl.at/scene/44eaf25f](https://superspl.at/scene/44eaf25f)
- mature and immature SILTs in normal colon: [https://superspl.at/scene/9a6e2550](https://superspl.at/scene/9a6e2550)

## ACCESS THE DATASET

### Data Availability

Dataset access details will be added here as data deposition is finalized. This section can be updated with direct links to image repositories, cloud buckets, portal records, and online viewers.

### Access Notes

Add access instructions here once the data hosting location and any access requirements are finalized.

## FILE ORGANIZATION

### Naming Conventions

Each file is expected to correspond to a specimen, sample, or derived analysis output. Naming conventions and release locations can be updated here as the dataset is finalized.

### File Types

| File Type | Description | Location |
| --- | --- | --- |
| `SampleID.ome.tif` | Stitched cyclic light-sheet image pyramid in OME-TIFF format | `TBD` |
| `SampleID-histology.ome.tif` | Histology image of matched or adjacent section, if available | `TBD` |
| `SampleID-mask.ome.tif` | Segmentation mask image | `TBD` |
| `SampleID-cells.csv` | Single-cell feature table including marker intensities and morphology | `TBD` |
| `SampleID-metadata.csv` | Sample- or region-level metadata | `TBD` |
| `SampleID-analysis.*` | Downstream analysis outputs and summaries | `TBD` |

## REPOSITORY LINKS

### Archive Links

- Zenodo archive: [https://doi.org/10.5281/zenodo.20170967](https://doi.org/10.5281/zenodo.20170967)
- CyANTs 3D registration workflow: [CyANTs README](CyANTs%203D%20registration/README.md)

## FILE LIST

### Current Repository Contents

- `README.md`: Project landing page and dataset information template.
- `CyANTs 3D registration/`: ANTsPy-based 3D registration workflow for large multichannel CyCIF / light-sheet Imaris volumes.
- `3D printed sample holders/`: 3D-printable circular spacers and glass stamps for mounting FFPE tissue in 6-well and 12-well glass-bottom plates.

## ADDITIONAL NOTES/COMMENTS

### Pending Updates

- Replace all `TBD` entries once the manuscript DOI, data portals, and archive records are available.
- Expand the file table with finalized naming conventions and storage locations before public release.
