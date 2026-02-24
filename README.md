# Project Kaiju #

### Introduction ###

`Kaiju` software includes the [Multiscale Atmosphere-Geospace Environment
(MAGE)](https://cgs.jhuapl.edu/Models/) model developed by the
[Center for Geospace Storms (CGS)](https://cgs.jhuapl.edu/) as well as
other scientific software for simulation of heliospheric environments
such as planetary magnetospheres and the solar wind. 

### Current version ###

Currently
supported applications include 
[MAGE](https://cgs.jhuapl.edu/Models/) (version 1.25)
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

⚠️ Note: The `master` branch (version 1.25) contains bleeding-edge code that is under
active development. If you are an inexperienced user, you might want
to stick with more stable [MAGE
0.75](https://github.com/JHUAPL/kaiju/tree/MAGE_0.75) and [MAGE
1.0](https://github.com/JHUAPL/kaiju/tree/MAGE_1.0) branches.

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

** GAMERA magnetosphere **: Sorathia, K. A., V. G. Merkin, E. V. Panov, B. Zhang,
J. G. Lyon, & J. Garretson, et al. (2020). Ballooning-interchange
instability in the near-Earth plasma sheet and auroral beads: Global
magnetospheric modeling at the limit of the MHD
approximation. Geophysical Research Letters, 47,
e2020GL088227. https://doi.org/10.1029/2020GL088227.

**REMIX:** Merkin, V. G., and J. G. Lyon (2010), Effects of the
low-latitude ionospheric boundary condition on the global
magnetosphere, J. Geophys. Res., 115, A10202.
https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010JA015461.

**RAIJU:** Team is working on the publication. 

**Dragon King:** Team is working on the publication.

**TIEGCM:** Qian, L., Burns, A.G., Emery, B.A., Foster, B., Lu, G.,
Maute, A., Richmond, A.D., Roble, R.G., Solomon, S.C. and Wang,
W. (2014). The NCAR TIE-GCM. In Modeling the Ionosphere–Thermosphere
System (eds J. Huba, R. Schunk and
G. Khazanov). https://doi.org/10.1002/9781118704417.ch7.

##### For the coupled MAGE model please cite the appropriate component models above and one of the following:

**MAGE 0.75:** Sorathia, K. A., Michael, A., Merkin,
V. G., Ohtani, S., Keesee, A. M., Sciola, A., et
al. (2023). Multiscale magnetosphere-ionosphere coupling during
stormtime: A case study of the dawnside current wedge. Journal of
Geophysical Research: Space Physics, 128,
e2023JA031594. https://doi.org/10.1029/2023JA031594.

**MAGE 1.0:** Pham, K. H., Zhang, B., Sorathia, K., Dang, T., Wang, W.,
Merkin, V., et al. (2022). Thermospheric density perturbations
produced by traveling atmospheric disturbances during August 2005
storm. Journal of Geophysical Research: Space Physics, 127,
e2021JA030071. https://doi.org/10.1029/2021JA030071.

**MAGE 1.25:** Team is working on the publication.


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
