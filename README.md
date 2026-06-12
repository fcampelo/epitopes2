# epitopes package
<!-- badges: start -->
  [![R-CMD-check](https://github.com/fcampelo/epitopes/workflows/R-CMD-check/badge.svg)](https://github.com/fcampelo/epitopes/actions)
  <!-- badges: end -->

Processing, Feature Extraction and modelling of epitope data from the Immune Epitope Database ([IEDB](http://iedb.org))

## Authors: 
- Felipe Campelo - [f.campelo@aston.ac.uk](mailto:f.campelo@aston.ac.uk), [fcampelo@gmail.com](mailto:fcampelo@gmail.com)

## Collaborators:
- Francisco Lobo - [franciscolobo@ufmg.br](mailto:franciscolobo@ufmg.br)

## Previous contributors:
- Jodie Ashford - [ashfojsm@aston.ac.uk](mailto:ashfojsm@aston.ac.uk)
- Joao Reis-Cunha - [jaumlrc@gmail.com](jaumlrc@gmail.com)

## Current version:
- Consolidation of linear B-cell epitopes from the full IEDB database and 
retrieval of protein data from NBCI and Uniprot.
- Organism- and taxon-based filtering of epitope data.
- Data splitting for machine learning based on protein dissimilarity.
- Extraction of observations for feature calculation using a sliding window representation.
- Calculation of sequence-based features.
- Read prediction output from other predictors (ABCPred, Bepipred 2.0, iBCE-EL, LBtope, SVMTrip)
- Fitting classification models to data with hold-out set-based assessment.

## Upcoming:
- Plotting functions

## Usage  
Check Vignette [Usage example](https://fcampelo.github.io/epitopes/articles/Usage-example.html) for a quick tutorial.

***

Please report any bugs or suggestions directly in the repository [Issues](https://github.com/fcampelo/epitopes/issues) page.
