# Project Kaiju #

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.16818620-blue)](https://doi.org/10.5281/zenodo.16818620)
[![Docs](https://img.shields.io/badge/Docs-readthedocs-blue)](https://kaiju-docs.readthedocs.io/en/latest/)

### Introduction ###

`Kaiju` software includes the [Multiscale Atmosphere-Geospace Environment
(MAGE)](https://cgs.jhuapl.edu/Models/) model developed by the
[Center for Geospace Storms (CGS)](https://cgs.jhuapl.edu/) as well as
other scientific software for simulation of heliospheric environments
such as planetary magnetospheres and the solar wind. 

### Current version ###

Currently
supported applications include 
[MAGE](https://cgs.jhuapl.edu/Models/) (version 1.5)
and [GAMERA-helio](https://cgs.jhuapl.edu/Models/gamera.php), i.e., the
geospace and inner heliosphere applications of the `kaiju`
software. `MAGE 1.25` includes [GAMERA](https://cgs.jhuapl.edu/Models/gamera.php) global magnetosphere
model, the `RAIJU` model of the inner magnetosphere, which is a full
rewrite of the [Rice Convection Model
(RCM)](https://cgs.jhuapl.edu/Models/rcm.php), the new energetic
precipitation model `Dragon King`, the [REMIX](https://cgs.jhuapl.edu/Models/remix.php) ionospheric electrodynamics
module and the
[TIEGCM](https://www.hao.ucar.edu/modeling/tgcm/tie.php) model of the
ionosphere-thermosphere.

⚠️ Note: The `master` branch (version 1.5) contains bleeding-edge code that is under
active development. If you are an inexperienced user, you might want
to stick with more stable [MAGE
1.25](https://github.com/JHUAPL/kaiju/tree/MAGE_1.25).

### NASA Community Coordinated Modeling Center (CCMC) availability ###

Users are also welcome to run previous versions of `MAGE` via the
NASA's [Community Coordinated Modeling
Center](https://ccmc.gsfc.nasa.gov/). `MAGE` versions
[0.75](https://ccmc.gsfc.nasa.gov/models/MAGE~0.75/) and
[1.0](https://ccmc.gsfc.nasa.gov/models/MAGE~1.0/) are available for
runs on request. `MAGE 0.75` includes `GAMERA`, `REMIX` and `RCM`. `MAGE
1.0` adds two-way coupling with `TIEGCM` to `MAGE 0.75`.

`GAMERA-helio` is also available for runs on request at the [CCMC](https://ccmc.gsfc.nasa.gov/models/GAMERA-Helio~1/).

### Documentation ###

Current documentation for the `kaiju` software is available via [Read The
Docs](https://kaiju-docs.readthedocs.io/en/latest/).

### Analysis ###

You are encouraged to use the [Kaipy](https://github.com/jhuapl/kaipy) package for analysis and
visualization of `Kaiju` simulations.

### Rules of the road ###

All users are strongly encouraged to contact the developers before
publication or presentation of `MAGE` or other `Kaiju` results. The
developers do their best to identify issues and bugs, but make no
guarantees. The developers are happy to help the users learn how to
apply the model correctly and interpret the results appropriately to
ensure their scientific robustness. Additional rules of the road can
be found in the [`Kaiju`
documentation](https://kaiju-docs.readthedocs.io/en/latest/roadrules.html).


### How to cite this work ###

We ask that the following papers be cited dependent on which
configuration the `Kaiju` software is being run.

#### For MAGE ####

##### For individual model components please cite the following:

**GAMERA MHD algorithm and tests:** Zhang, B., Sorathia, K.A., Lyon, J.G., Merkin, V.G., Garretson,
J.S. and Wiltberger, M., 2019. GAMERA: A three-dimensional
finite-volume MHD solver for non-orthogonal curvilinear
geometries. The Astrophysical Journal Supplement Series, 244(1),
p.20. https://iopscience.iop.org/article/10.3847/1538-4365/ab3a4c/meta.

**GAMERA magnetosphere**: Sorathia, K. A., V. G. Merkin, E. V. Panov, B. Zhang,
J. G. Lyon, & J. Garretson, et al. (2020). Ballooning-interchange
instability in the near-Earth plasma sheet and auroral beads: Global
magnetospheric modeling at the limit of the MHD
approximation. Geophysical Research Letters, 47,
e2020GL088227. https://doi.org/10.1029/2020GL088227.

**REMIX:** Merkin, V. G., and J. G. Lyon (2010), Effects of the
low-latitude ionospheric boundary condition on the global
magnetosphere, J. Geophys. Res., 115, A10202.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JA015461.

**RAIJU:** Sciola, A., Sorathia, K., Merkin, V. G., Shumko, M.,
Martinis, C. R., & Garretson, J. (2025). Connecting stormtime SAR arc
variability to mesoscale ring current-plasmasphere evolution via
global geospace modeling. AGU Fall Meeting 2025, SM51B-114.
https://ui.adsabs.harvard.edu/abs/2025AGUFMSM51B.114S/abstract.

**Dragon King:** Lin, D., Wang, W., Merkin, V. G., Wiltberger, M. J.,
Bao, S., Sorathia, K., et al. (2025). Characterizing auroral
precipitation and energy deposition with the Dragon King model in
MAGE. AGU Fall Meeting 2025, SA51C-1132.
https://ui.adsabs.harvard.edu/abs/2025AGUFMSA51C1132L/abstract.

Bao, S., Lin, D., Sorathia, K., Wang, D., Toffoletto, F., Merkin,
V. G., et al. (2022). Wave-dependent ring current electron loss and
magnetospheric precipitation during a geospace storm. AGU Fall
Meeting 2022, SM22F-1988.
https://ui.adsabs.harvard.edu/abs/2022AGUFMSM22F1988B/abstract.

**TIEGCM:** Wu, H., Wang, W., Pham, K. H., Lin, D., Rao, N.,
Wiltberger, M. J., et al. (2025). The NCAR-TIEGCM Version 3.0.
Journal of Geophysical Research: Space Physics, 130(9),
e2025JA034219. https://doi.org/10.1029/2025JA034219.

##### For the coupled MAGE model please cite the appropriate component models above and one of the following:

**MAGE software, all versions:** https://doi.org/10.5281/zenodo.16818620.

**MAGE 0.75:** Sorathia, K. A., Michael, A., Merkin,
V. G., Ohtani, S., Keesee, A. M., Sciola, A., et
al. (2023). Multiscale magnetosphere-ionosphere coupling during
stormtime: A case study of the dawnside current wedge. Journal of
Geophysical Research: Space Physics, 128,
e2023JA031594. https://doi.org/10.1029/2023JA031594.
Software: MAGE_0.75.5, https://doi.org/10.5281/zenodo.16818621.

**MAGE 1.0:** Pham, K. H., Zhang, B., Sorathia, K., Dang, T., Wang, W.,
Merkin, V., et al. (2022). Thermospheric density perturbations
produced by traveling atmospheric disturbances during August 2005
storm. Journal of Geophysical Research: Space Physics, 127,
e2021JA030071. https://doi.org/10.1029/2021JA030071.
Software: MAGE_1.0.5, https://doi.org/10.5281/zenodo.16818682.

**MAGE 1.25:** Team is working on the publication.
Software: MAGE_1.25.1, https://doi.org/10.5281/zenodo.17915213
(MAGE_1.25.0, https://doi.org/10.5281/zenodo.17914755).

**MAGE 1.5:** Team is working on the publication.

#### For GAMERA-helio ####

Provornikova, E., Merkin, V.G., Vourlidas, A., Malanushenko, A.,
Gibson, S.E., Winter, E. and Arge, C.N., 2024. MHD Modeling of a
Geoeffective Interplanetary Coronal Mass Ejection with the Magnetic
Topology Informed by In Situ Observations. The Astrophysical Journal,
977(1), p.106. https://iopscience.iop.org/article/10.3847/1538-4357/ad83b1/meta.

Merkin, V. G., J. G. Lyon, D. Lario, C. N. Arge, and C. J. Henney
(2016), Time-dependent magnetohydrodynamic simulations of the inner
heliosphere, J. Geophys. Res. Space Physics, 121, 2866–2890.
https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2015JA022200.


### License ###

`Kaiju` is distributed under the [BSD 3-Clause
license](LICENSE.md).


### Contribution guidelines ###

All contributions should be made by forking this repository and
  submitting a pull request.


### Who do I talk to? ###

Feel free to send us a message via the feedback form [here](https://cgs.jhuapl.edu/feedback/).
