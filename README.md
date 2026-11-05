# Dysregulated B Cell Education in Down Syndrome
Scripts used for the computational analysis in the Medrxiv preprint "Dysregulated B Cell Education Drives Autoimmunity and Immunodeficiency in Down Syndrome" (Ramba et al)

# Workflow Overview and Dependencies
Analysis is split into spatial transcriptomic analysis and single cell RNA/BCR sequencing analysis. Each section requires its own dependencies (available for installation via yml conda files for each section). 

The code is organized into several R notebooks, covering different steps of analysis, plus necessary additional files and shell scripts. 

Spatial Transcriptomics:
| SCRIPT | NOTEBOOK LINK | DESCRIPTION  |
|------------------|---------------------------|---------------------------|
| YML      | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/spatial_transcriptomics/spatial_environment.yml)               | YML file to recreate conda environment for spatial analysis |
| Gene expression analysis - Part 1 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/spatial_transcriptomics/spatial_GEX_analysis/Processing_1.R)| Processing of spatial h5 files |
| Gene expression analysis - Part 2 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/spatial_transcriptomics/spatial_GEX_analysis/GEX_analysis_2.R)| Analysis of processed files: Figure 2, Figure 3A-E, Supplemental Figure 2, 3 |
| Pseudobulking analysis - Part 1 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/spatial_transcriptomics/pseudobulking/Metadata.R)| Generate metadata file for pseudobulk analysis |
| Pseudobulking analysis - Part 2 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/spatial_transcriptomics/pseudobulking/Pseudobulking_input.R)| Generate input files for pseudobulk analysis |
| Pseudobulking analysis - Part 3 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/spatial_transcriptomics/pseudobulking/Pseudobulking_group+count.R)| Pseudobulk analysis: Figure 3F |

Single Cell RNA-sequencing/BCR-sequencing:
| SCRIPT | NOTEBOOK LINK | DESCRIPTION  |
|------------------|---------------------------|---------------------------|
| YML      | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/complete_env.yml) | YML file to recreate conda environment for scRNAseq+BCRseq analysis (NOTE: need to install additional packages, see README file within folder) |
| Sample Data | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/sample_data.csv)| Necessary CSV file for running analysis |
| Gene expression analysis - Part 1| [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/gex_analysis/B_GEX_analysis.R) | Processing of matrix files & integration |
| Gene expression analysis - Part 2 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/gex_analysis/B_GEX_analysis_2.R) | Analysis |
| Gene expression analysis - Part 3 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/gex_analysis/B_GEX_analysis_3.R) | Analysis of processed files: Figure 4B-D, Supplemental Figure 4 |
| Gene expression analysis - Part 4 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/gex_analysis/ucell_analysis.R) | UCell analysis of processed files: Figure 4F-G, Supplemental Figure 4|
| BCR analysis - Part 1 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/bcr_analysis/igblast.sh) | IgBlast analysis |
| BCR analysis - Part 2 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/bcr_analysis/tsv_merge.R) | Merges individual files into single file for downstream analysis|
| BCR analysis - Part 3 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/bcr_analysis/immcantation.R) | Immcantation analysis: Supplemental Figure 5 |
| BCR analysis - Part 4 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/bcr_analysis/immcantation_2.R) | Immcantation analysis: Figure 5D |
| BCR analysis - Part 5 | [link](https://github.com/meredithramba/DownSyndrome_Bcells/blob/main/scRNAseq%2BBCRseq/bcr_analysis/immcantation_3.R) | Immcantation analysis: Figure 4E, 5A-C, E, Supplemental Figure 4, 5|
