# SamsaRaLight

**SamsaRaLight** is an R package for simulating light transmission in forest stands based on the ray-tracing model SamsaraLight ([Courbaud et al. 2003](https://www.sciencedirect.com/science/article/pii/S016819230200254X?via%3Dihub)).

📊 It allows users to estimate:

-   Light intercepted by individual trees
-   Light reaching the forest floor
-   Light at virtual ground sensors


✨  With many features: 

-   Estimation of direct and diffuse radiation
-   Flexible crown geometry (symmetric and asymmetric)
-   Multiple transmission models
-   Light competition indicators
-   Virtual light sensors
-   2D visualization tools
-   Parallel simulations


🌐 Model documentation and tutorials are available here: <https://natheob.github.io/SamsaRaLight/>


## Installation

You can install the development version of SamsaRaLight from [GitHub](https://github.com/NatheoB/SamsaRaLight/) with:

``` r
install.packages("devtools")
devtools::install_github("NatheoB/SamsaRaLight")
```

## Citation

If you use the SamsaraLight in your research, please cite :

- The original publication of the SamsaraLight model: Courbaud, B., de Coligny, F. & Cordonnier, T. Simulating radiation distribution in a heterogeneous Norway spruce forest on a slope. Agricultural and Forest Meteorology 116, 1–18 (2003).
[Courbaud et al. 2003](https://www.sciencedirect.com/science/article/pii/S016819230200254X?via%3Dihub).
- This package using `r citation("SamsaRaLight")`: [Beauchamp et al. 2026](https://doi.org/10.32614/CRAN.package.SamsaRaLight).


## Contributing

Contributions, bug reports, and feature requests are welcome.

Please [open an issue](https://github.com/NatheoB/SamsaRaLight/issues) or [pull request](https://github.com/NatheoB/SamsaRaLight/pulls) on the github project.

## License

This package is released under the GPL-3 License. See LICENSE for details.
