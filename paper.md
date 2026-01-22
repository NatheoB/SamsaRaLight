# NA

## Summary

The *SamsaRaLight* R package (Figure 1) provides a simplified and
efficient R interface to the SamsaraLight ray-tracing model
\[@courbaud_simulating_2003\], enabling explicit simulation of light
competition in forest stands. It allows users to estimate the energy
absorbed by individual trees, the light transmitted to the forest floor,
or the radiation reaching user-defined locations via virtual sensors.
Forest stands are represented in 3D, with individual crowns modelled
using symmetric or asymmetric geometric shapes and stand attributes
including size, slope, and aspect. Direct and diffuse radiation is
discretised according to location and cast through the canopy, allowing
light attenuation across successive crowns to be explicitly simulated.
Outputs include light interception, transmission, and competition
variables at tree, ground, or sensor levels, supporting studies of tree
growth, mortality, and regeneration dynamics within reproducible and
optimised R workflows.

![Global aim of the SamsaRaLight R
package](paper/figures/fig_global.jpg)

Global aim of the SamsaRaLight R package

## Statement of need

### Light competition and the need for ray-tracing approaches

Light is a key resource driving the dynamics of saplings and trees
\[@binkley_light_2013\], inducing strong plasticity in crown dimensions
under competitive pressure \[@touzot_shade_2025\]. Numerous forest
radiative transfer models have therefore been developed to represent
canopy structure and light absorption processes \[@ligot_forest_2014\].
A defining feature of light competition is its strong asymmetry,
resulting from the attenuation of light along its path through
successive canopy layers \[@weiner_asymmetric_1990;
@schwinning_mechanisms_1998\]. Accurately representing this process
requires ray-tracing approaches that explicitly account for light path
length across overlapping tree crowns.

The SamsaraLight ray-tracing model, developed by Courbaud et
al. \[-@courbaud_simulating_2003\], is one of the most widely used
models for this purpose in forest ecology. It has been applied to study
the effects of light interception on tree growth and mortality dynamics
\[@beauchamp_light_2025\], as well as to quantify light transmittance to
the forest floor and its influence on regeneration processes
\[@ligot_managing_2014; @ligot_tree_2016\].

### Integration within forest simulators and limits of existing implementations

Explicitly coupling light competition with forest simulators is
essential to analyse how silvicultural strategies influence stand
structure, dynamics, and regeneration \[@lafond_uneven-aged_2014\]. For
this reason, SamsaraLight has been integrated as a core module in
numerous forest simulators developed within the CAPSIS platform
\[@dufour-kowalski_capsis_2012\], including Samsara2
\[@courbaud_applying_2015\], PDG-ARENA \[@rouet_pdg-arena_2025\],
HETEROFOR 1.0 \[@jonard_heterofor_2020; @de_wergifosse_heterofor_2020\],
and RReShar \[@barrere_oak_2024\].

However, the tight coupling of SamsaraLight to CAPSIS, a Java-based
framework, limits its independent use and integration into common
statistical workflows. To address this, Fortin
\[-@fortin_executing_2020\] introduced the R packages J4R and RCapsis,
enabling CAPSIS models to be accessed from R. Despite this advance, the
approach remains technically demanding, requiring interaction with
external software, intermediate files, and multi-language pipelines,
which restricts its accessibility for many users.

### Simplifying and broadening access through an R implementation

Many researchers addressing forest dynamics and regeneration questions
increasingly rely on R for model fitting, simulation analysis, and
reproducible workflows, often outside the CAPSIS ecosystem. While some
forest models have been fully reimplemented in R (e.g. SurEau-Ecos
\[@ruffault_sureau-ecos_2022\]), ray-tracing approaches remain
computationally intensive and difficult to integrate without complex
software dependencies. A native R implementation of SamsaraLight
simplifies the explicit representation of light interception and
transmission, enabling its direct use in studies of tree growth,
mortality, and regeneration.

Beyond simplifying existing workflows, the SamsaRaLight R package
facilitates links with emerging technologies and research fields. In
particular, the increasing availability of LiDAR data in forest ecology
\[@lines_shape_2022\] creates opportunities to couple ray-tracing light
models with LiDAR-based crown reconstruction tools such as ITSMe
\[@terryn_analysing_2023\] and LidaRtRee \[@monnet_lidartree_2025\].
More broadly, providing an accessible R implementation may foster
applications beyond forestry, including agroforestry, where tree shading
strongly influences crop production \[@dupraz_hi-safe_2019\].

## Key features

### Native R workflow with pre- and post-processing tools

The *SamsaRaLight* R package provides a fully native R implementation of
the SamsaraLight ray-tracing model. Forest stands and radiation inputs
are defined using standard R formats, with dedicated functions to check
and construct virtual stands from tree inventories and to derive monthly
radiation data based on stand location. Monthly global radiation and
diffuse-to-global ratios are retrieved from the European PVGIS database
\[@huld_new_2012\]. Simulations are executed with default settings
suitable for most applications, while allowing advanced users to adjust
ray-tracing parameters and internal model behaviour. Standard R methods
for S3 objects ([`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)) enable
inspection and visualisation of both input stands and output light
variables, with graphical outputs generated using *ggplot2*
\[@wickham_ggplot2_2016\].

### High computational performance through C++ implementation and parallelisation

To ensure efficient computation, the core of the SamsaraLight model is
implemented in C++ and interfaced with R using *Rcpp*
\[@eddelbuettel_rcpp_2025\]. This design leverages the computational
performance and memory management of C++ while preserving the
flexibility of R-based workflows. Ray-tracing calculations can be
executed in parallel using shared-memory parallelisation with OpenMP,
with fine control over the number of cores and the ability to activate
or deactivate parallel execution depending on computational constraints.

### Documentation and learning resources

The *SamsaRaLight* package is accompanied by an online documentation
website generated with *pkgdown* \[@wickham_pkgdown_2025\], available at
<https://natheob.github.io/SamsaRaLight/>. The documentation provides a
detailed description of the underlying model and a set of progressive
vignettes designed to support learning and practical use. These
tutorials are based on real forest datasets distributed with the package
and guide users through typical workflows and parameter choices.

## Acknowledgments

The SamsaraLight model has been initially designed by Benoit Courbaud
and implemented in Capsis by François de Coligny. SamsaraLight has been
next updated with new functionality or upgraded for faster computation
time by various persons across years.
