# Interspecific-hybridisation-provides-a-low-risk-option for increasing genetic diversity of reef-building corals
>This repository contains the code and data required to reproduce analyses conducted by Annika Lamb (AIMS) for the manuscript "Interspecific hybridisation provides a low risk option for increasing genetic diversity of reef-building corals".
>
>Data is avaiable in .csv, .fasta, and .txt formats and code is written in R.
>
# Data
>Repository contains five phenotypic datasets and one ITS2 metabarcoding dataset
>
## Fertilisation data
>HybridFieldProject_FertilisationData.csv
>
>Fertilisation success of purebred and hybrid crosses.
>
>Dataset comprised of:
- Mother: dam ID
+ Cross: offspring Group as Acropora sarmentosa purebred (SS), Acropora florida purebred (FF), FS hybrid, or SF hybrid
* Sample: replicate count
- SampleID: unique sample ID
+ Fertilised: number of the 100 eggs that had divided
* Unfertilised: number of the 100 eggs that had not divided
- Proportion: proportion of divided eggs
+ Percentage: percentage of divided eggs
## Survival data
### KM dataset
>HybridFieldExperiment_Survival_KMDataset.csv
>
>Survival of the coral juveniles over time in format appropriate for Kaplan-Meier Estimations.
>
>Dataset comprised of:
- TileID: Unique Tile ID
+ Cross: Offspring Group as Acropora sarmentosa purebred (SS), Acropora florida purebred (FF), FS hybrid, or SF hybrid
* Rod: number of the rod the coral was assigned to
- Cassette: number of the cassette the coral was assigned to
+ Frame: number of the frame the coral was assigned to
* RecruitID: recruit number on tile
- UniqueID: unique recruit ID
+ Time: time of sampling in months
* Days: time of sampling in days
- Date: date of sampling
+ Deployment: date of deployment
* Survival: survival where '1' is dead and '2' is alive
### BGLMM dataset
>HybridFieldExperiment_Survival_BGLMMDataset.csv
>
>Survival of the coral juveniles over time in format appropriate for Bayes Generalised Linear Mixed Effects Modelling.
>
>Dataset comprised of:
- TileID: Unique Tile ID
+ Cross: Offspring Group as Acropora sarmentosa purebred (SS), Acropora florida purebred (FF), FS hybrid, or SF hybrid
* Rod: number of the rod the coral was assigned to
- Cassette: number of the cassette the coral was assigned to
+ Frame: number of the frame the coral was assigned to
* RecruitID: recruit number on tile
- UniqueID: unique recruit ID
+ Location: position on tile
* TimePoint: time of sampling in months
- Date: date of sampling
+ Survival: survival where '1' is alive and '0' is dead
* Month: month of sampling
## Size data
>HybridFieldExperiment_GrowthDataset.csv
>
>Size of the coral juveniles over time.
>
>Dataset comprised of:
- TileID: Unique Tile ID
+ Cross: Offspring Group as Acropora sarmentosa purebred (SS), Acropora florida purebred (FF), FS hybrid, or SF hybrid
* Rod: number of the rod the coral was assigned to
- Cassette: number of the cassette the coral was assigned to
+ Frame: number of the frame the coral was assigned to
* RecruitID: recruit number on tile
- UniqueID: unique recruit ID
+ Location: position on tile
* TimePoint: time of sampling in months
- Date: date of sampling
+Size: surface area of the recruit in mm^2
## Colour data
>HybridFieldExperiment_ColourDataset.csv
>
>Colour of the coral juveniles over time.
>
>Dataset comprised of:
- TileID: unique Tile ID
+ Cross: offspring Group as Acropora sarmentosa purebred (SS), Acropora florida purebred (FF), FS hybrid, or SF hybrid
* Rod: number of the rod the coral was assigned to
- Cassette: number of the cassette the coral was assigned to
+ Frame: number of the frame the coral was assigned to
* RecruitID: recruit number on tile
- UniqueID: unique recruit ID
+ Location: position on tile
* TimePoint: time of sampling in months
- Date: date of sampling
+ Survival: survival where '1' is alive and '0' is dead
* Colour: calibrated colour value relative to the Coral Watch Colour Chart
- RelativeGrey: grey value relative to starting colour
+ Month: month of sampling
## ITS2 metabarcoding
> FieldITS2.profiles.absolute.abund_and_meta.txt: sequence/profile count data
> 
> FieldITS2.seqs.fasta: sequence data
> 
> metadata.csv:
- sample unique ID (sample_uid)
+ sample name (sample_name)
* Preservation method (ethanol or liquid nitrogen)
- Offspring group (Species)
+ Frame number (Frame)
* Tile ID (Tile)
- Tile number (TileNumber)
+ Sample type (ControlSample: control or sample)
* Recruit size in mm^2 (Size)

# Code
>Repository contains two scripts for analysing and graphing phenotypic and ITS2 data of hybrid and purebred corals
>
## Phenotypic data analysis
>R code for analysing hybrid and purebred phenotypic data.
>
>Lamb_etal_BioOpen_Code.R
>
>Requires HybridFieldExperiment_ColourDataset.csv, HybridFieldExperiment_GrowthDataset.csv, HybridFieldExperiment_Survival_BGLMMDataset.csv, HybridFieldExperiment_Survival_KMDataset.csv, and HybridFieldProject_FertilisationData.csv
>
## ITS2 metabarcoding
>R code for analysing ITS2 metabarcoding data.
>
>Lamb_etal_BioOpen_ITS2metabarcoding.R
>
>Requires FieldITS2.profiles.absolute.abund_and_meta.txt, FieldITS2.seqs.fasta, and metadata.csv
