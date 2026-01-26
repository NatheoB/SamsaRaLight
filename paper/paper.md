---
title: "SamsaRaLight: An R package for estimation of tree light interception using ray-tracing"
tags:
  - R
  - forest modelling
  - ray tracing
  - light competition
  - ecology
authors:
  - name: Nathéo Beauchamp
    affiliation: 1
    orcid: 0009-0007-9103-5194
  - name: François de Coligny
    affiliation: 2
    orcid: 0000-0002-8538-3009
  - name: Benoit Courbaud
    affiliation: 3
    orcid: 0000-0002-3050-9559
  - name: Maxime Jaunatre
    affiliation: 3
    orcid: 0009-0002-2816-1677
  - name: Gauthier Ligot
    affiliation: 1
    orcid: 0000-0002-5508-4358
affiliations:
  - name: University of Liege, Gembloux Agro-Bio Tech, Unité de Gestion des Ressources forestières et des Milieux naturels, B-5030 Gembloux, Belgium
    index: 1
  - name: botAnique et Modélisation de l’Architecture des Plantes et des végétations AMAP, Montferrier-sur-Lez, France
    index: 2
  - name: University of Grenoble Alpes, LESSEM, INRAE, Grenoble, France
    index: 3
bibliography: paper.bib
---

## Summary

The *SamsaRaLight* R package (Figure 1) provides a simplified and efficient R interface to the SamsaraLight ray-tracing model [@courbaud_simulating_2003], enabling explicit simulation of light competition in forest stands. It allows users to estimate the energy absorbed by individual trees, the light transmitted to the forest floor, and optionally the radiation reaching user-defined locations via virtual sensors. Forest stands are represented in 3D, with individual crowns modelled using symmetric or asymmetric geometric shapes and stand attributes including size, slope, and aspect. Direct and diffuse radiation is discretised according to stand location and cast through the canopy, allowing light attenuation across successive crowns to be explicitly simulated. Outputs include light interception, transmission, and competition variables at tree, ground, and sensor levels, supporting studies of tree growth, mortality, and regeneration dynamics within reproducible and optimised R workflows.

![Global aim of the SamsaRaLight R package. The example pipeline is executed on the stored dataset `SamsaRaLight::data_cloture20`](figures/fig_global.jpg)

## Statement of need

Light is a key resource driving the dynamics of saplings and trees [@binkley_light_2013], inducing strong plasticity in crown dimensions under competitive pressure [@touzot_shade_2025]. Numerous forest radiative transfer models have therefore been developed to represent canopy structure and light absorption processes [@ligot_forest_2014]. A key feature of light competition is its strong asymmetry, resulting from the attenuation of light along its path through successive canopy layers [@weiner_asymmetric_1990; @schwinning_mechanisms_1998]. Accurately representing this process requires ray-tracing approaches that explicitly account for light path length across overlapping tree crowns.

The SamsaraLight ray-tracing model, developed by Courbaud et al. [-@courbaud_simulating_2003], is one of the most widely used models for this purpose in forest ecology. It has been applied to study the effects of light interception on tree growth and mortality dynamics [@beauchamp_light_2025], as well as to quantify light transmittance to the forest floor and its influence on regeneration processes [@ligot_managing_2014; @ligot_tree_2016]. Explicitly coupling light competition with forest simulators is essential to study forest stand evolution over time and how silvicultural strategies influence stand structure, dynamics, and regeneration [@lafond_uneven-aged_2014, @barrere_oak_2024]. For this reason, SamsaraLight has been integrated as a core module in numerous forest simulators including Samsara2 [@courbaud_applying_2015], PDG-ARENA [@rouet_pdg-arena_2025], HETEROFOR 1.0 [@jonard_heterofor_2020; @de_wergifosse_heterofor_2020], or RReShar [@barrere_oak_2024].

### 

## State of the field

The SamsaraLight ray-tracing model was originally developed within the CAPSIS platform [@dufour-kowalski_capsis_2012], a Java-based framework simplifying the develoment of forest simulators. However, the tight coupling of SamsaraLight to CAPSIS limits its independent use and integration into common statistical workflows. To address this, Fortin [-@fortin_executing_2020] introduced the R packages J4R and RCapsis, enabling CAPSIS models to be accessed from R. Despite this advance, the approach remains technically demanding, requiring interaction with external software, intermediate files, and multi-language pipelines, which restricts accessibility for many users.

Many researchers working on forest dynamics and regeneration increasingly rely on R for model fitting, simulation analysis, and reproducible workflows, often outside the CAPSIS ecosystem. While some forest models have been fully reimplemented in R (*e.g.*, SurEau-Ecos [@ruffault_sureau-ecos_2022]), ray-tracing approaches remain computationally intensive and difficult to use efficiently without high-performance implementations and streamlined interfaces. By combining a user-friendly R interface with an efficient C++ backend, the *SamsaRaLight* R package addresses these limitations, simplifying the explicit representation of light interception and transmission and enabling direct use in studies of tree growth, mortality, and regeneration.

## Software design

### Native R workflow with pre- and post-processing tools

The *SamsaRaLight* R package provides a fully native R implementation of the SamsaraLight ray-tracing model. Forest stands and radiation inputs are defined using standard R formats, with dedicated functions to check and construct virtual stands from tree inventories and to derive monthly radiation data based on stand location. Monthly global radiation and diffuse-to-global ratios are retrieved from the European PVGIS database [@huld_new_2012]. Simulations are executed with default settings suitable for most applications, while allowing advanced users to adjust ray-tracing parameters and internal model behaviour. Standard R methods for S3 objects (`print()`, `summary()`, and `plot()`) enable inspection and visualisation of both input stands and output light variables, with graphical outputs generated using *ggplot2* [@wickham_ggplot2_2016].

### High computational performance through C++ implementation and parallelisation

To ensure efficient computation, the core of the SamsaraLight model is implemented in C++ and interfaced with R using *Rcpp* [@eddelbuettel_rcpp_2025]. This design leverages the computational performance and memory management of C++ while preserving the flexibility of R-based workflows. Ray-tracing calculations can be executed in parallel using shared-memory parallelisation with OpenMP, with fine control over the number of cores and the ability to activate or deactivate parallel execution depending on computational constraints.

### Documentation and learning resources

The *SamsaRaLight* package is accompanied by an online documentation website generated with *pkgdown* [@wickham_pkgdown_2025], available at <https://natheob.github.io/SamsaRaLight/>. The documentation provides a detailed description of the underlying model and a set of progressive vignettes designed to support learning and practical use. These tutorials are based on real forest datasets distributed with the package and guide users through typical workflows and parameter choices.

## Research impact statement

The *SamsaRaLight* R package has multiple utilities for researchers, students, and field foresters, enabling broader and more user-friendly access to ray-tracing models. Its applications include: (1) simplifying workflows for forest research involving light competition; (2) serving as a teaching tool through progressive tutorials that guide students and foresters in understanding and applying ray-tracing models; and (3) providing accessible visualization tools for estimating light distribution in forest stands directly from tree inventories, with future plans to extend usability via a Shiny application for non-programmer users [@shiny_2025].

The *SamsaRaLight* R package could facilitate links with emerging technologies and research fields. In particular, the increasing availability of LiDAR data in forest ecology [@lines_shape_2022] creates opportunities to couple ray-tracing light models with LiDAR-based crown reconstruction tools such as ITSMe [@terryn_analysing_2023] and LidaRtRee [@monnet_lidartree_2025]. More broadly, providing an accessible R implementation may foster applications beyond forestry, including agroforestry, where the SamsaraLight model has already been used for estimating tree shading influences on crop production [@dupraz_hi-safe_2019].

## AI usage disclosure

Generative AI tools were used in a limited support role. The core C++ implementation of the SamsaRaLight model and the R interface were fully developed by the authors before the widespread availability and use of generative AI tools. ChatGPT was later used to help refine function documentation, tests, and user-facing interfaces for inspection. The initial manuscript and online documentation were first written by the authors without AI assistance, and ChatGPT and DeepL Write were subsequently used only to improve clarity, conciseness, and reduce redundancy in the text. All AI-assisted content was reviewed, edited, and validated by the authors, who take full responsibility for the software and the paper.

## Acknowledgments

The SamsaraLight model was initially designed by Benoit Courbaud (University of Grenoble Alpes, LESSEM, France) and François de Coligny (University of Montpellier, AMAP, France), and subsequently developed into a library within the CAPSIS Java platform with contributions from Nicolas Donès (INRAE PIAF, Clermont-Ferrand, France), Gauthier Ligot (University of Liège, Gembloux, Belgium), Mathieu Jonard, and Frédéric André (UCL, Belgium), who also added new functionality and improved computational performance over the years. No dedicated funding was used for the development of the *SamsaRaLight* R package, which was conducted independently as a personal side project. We thank the University of Grenoble Alpes (France) and the University of Liège (Belgium) for providing a supportive research environment while this work was carried out alongside the main research activities of Nathéo Beauchamp.
