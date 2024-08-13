# Interspecific-hybridisation-provides-a-low-risk-option for increasing genetic diversity of reef-building corals
>This repository contains the code and data required to reproduce analyses conducted by Annika Lamb (AIMS) for the manuscript "Interspecific-hybridisation-provides-a-low-risk-option for increasing genetic diversity of reef-building corals".
>
>Data is avaiable in .csv, .fasta, and .txt formats and code is written in R.
>
# Data
>Repository contains four phenotypic datasets, one ITS2 metabarcoding dataset, and one temperature profile dataset
>
## Fertilisation data
>
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
>
## Size data
>HeatWaveExperiment_Results_Singles_Growth.csv
>
>Size of the coral juveniles over time.
>
>Dataset comprised of:
- TileID: Unique Tile ID
+ Cross: hybrid (TL) or purebred (LL or TT) experimental group
* Batch: fertilisation cross (TT1, TT2, TL1, TL2, LL1, or LL2)
- Tank: holding tank
+ Recruit: number of recruit on tile
* RunningRecruit: unique recruit number
- Size: mm^2
+ Timepoint: sampling timepoint
* Day: sampling day
- DHW: number of accumulated degree heating weeks of thermal stress
+ Treatment: temperature treatment as ambient or elevated
>
## Colour data
>HeatWaveExperiment_Results_Singles_Colour.csv
>
>Colour of the coral juveniles over time.
>
>Dataset comprised of:
- TileID: Unique Tile ID
+ Cross: hybrid (TL) or purebred (LL or TT) experimental group
* Batch: fertilisation cross (TT1, TT2, TL1, TL2, LL1, or LL2)
- Tank: holding tank
+ Recruit: number of recruit on tile
* RunningRecruit: unique recruit number
- Colour: recruit colour score
+ Timepoint: sampling timepoint
* Day: sampling day
- DHW: number of accumulated degree heating weeks of thermal stress
+ Treatment: temperature treatment as ambient or elevated
>
## ITS2 metabarcoding
> HeatWaveITS2.seqs.absolute.abund_and_meta.txt: sequence/profile count data
> 
> HeatWaveITS2.seqs.fasta: sequence data
> 
> metadata.csv:
- sample unique ID (sample_uid)
+ sample name (sample_name)
* cross (species)
- holding tank number (Tank)
+ tile ID (Tile)
* temperature treatment (Treatment: ambient or elevated)
- sample type (ControlSample: control or sample).

## Temperature profiling
> Temperature profile of the experiment
> Dataset comprised of:
- RampingDate: date of experiment
+ DaviesDate: date of temperature recording at Davies Reef by the Davies Reef Weather Station (Australian Institute of Marine Science)
* Day: day of the experiment
- ElevatedTemperature: temperature in the elevated treatment
+ AmbientTemperature: temperature in the ambient treatment
* RampTemp: temperature in elevated treatment post ramping to ambient temperature
- DegreesAboveAmbientTreatment: number of degrees above ambient in the elevated treatment
+ DegreesAboveDaviesFebruaryAverage: degrees above the maximum monthly mean at Davies Reef

# Code
>Repository contains three scripts for analysing and graphing phenotypic and ITS2 data of hybrid and purebred corals
>
## Phenotypic data analysis
>R code for analysing hybrid and purebred phenotypic data.
>
>Lamb_etal_PhenotypicAnalyses.R
>
>Requires HeatWaveExperiment_Results_Singles_Colour.csv, HeatWaveExperiment_Results_Singles_Growth.csv, HeatWaveExperiment_Results_Singles_Survival.csv, and HybridPAM_PrelimResults_RecruitCentres_T4_Long.csv
>
## ITS2 metabarcoding
>R code for analysing ITS2 metabarcoding data.
>
>Lamb_etal-ITS2_sequence_analysis.R
>
>Requires HeatWaveITS2.profiles.absolute.abund_and_meta.txt, HeatWaveITS2.seqs.fasta, and metadata.csv
>
## Fertilisation data analyses
>R code for analysing fertilisation data
>
>Lamb_etal_TemperatureProfileGraphing.R
>
>Requires FinalTemperatureProfile_DHWsIncluded.csv
