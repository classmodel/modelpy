[![CLASS4GL Logo](https://class4gl.eu/wp-content/uploads/2019/01/cropped-class4gl_small-1.png)](https://class4gl.eu)

_CLASS4GL_ (Chemistry Land-surface Atmosphere Soil Slab model for Global Studies) is a fast and easy interface to investigate the dynamics of the atmospheric boundary layer from weather balloons worldwide. General info and tutorials for using CLASS4GL are available at https://class4gl.eu, and video clips about the atmospheric boundary layer physics can be found on the [website of the original CLASS model](classmodel.github.io/).

# Features
  - _Mine_ appropriate observations from global radio soundings, satellite data, reanalysis and climate models
  - _Automise_ mass parallel simulations of the atmospheric boundary layer and global sensitivity experiments
  - _Foster_ a better understanding of land-atmosphere interactions and the drivers of extreme weather globally
  - _Share_ your data, experiments, and code developments with the research community

# Method

### Description

The framework CLASS4GL is designed to facilitate the investigation of the atmospheric boundary layer evolution in response to different land and atmospheric conditions observed around the world. The core of the platform is the model CLASS that is used to simulate the evolution of the atmospheric boundary layer. Instruction video about the boundary layer processes and how they are considered in the CLASS model can be found as on the [CLASS model website](https://classmodel.github.io/). Observational data from balloons, satellites and reanalysis, are used to constrain and initialize the model. CLASS4GL uses 2 million global balloon soundings from the integrated global radio sounding archive and satellite data from the last 40 years.

### Components

  - A global data module that employs balloon soundings, satellite imagery and reanalysis data
  - An interface to easily perform multiple simulations of the atmospheric boundary layer in parallel, and multiple batches of global sensitivity experiments
  - Tools for Pre-and post-processing the data pool of input data and experiments.
  - A GUI data explorer

The tool is under continuous development, and it can downloaded and installed as described in the tutorials on class4gl.eu/#getstarted.

In case you experience a problem or a bug, please don’t hesitate to contact us class4gl.eu/#contact. You an also open an issue on the github page (https://github.com/hendrikwout/class4gl/issues) . Any feedback will be highly appreciated.

### Data sources

CLASS4GL employs the balloon soundings from the Integrated Global Radiosonde Archive (IGRA) to initialize and validate the CLASS model. The sounding data is supplemented with ancillary data to further constrain the model. Therefore, a default set of gridded global datasets from satellite imagery, reanalysis and and surveys have been used that span a period of 1981–2015. An complete overview of the datasets can be found in the table. However, the default set can be replaced by alternative datasets as long as they are provided in netCDF format.

[Schematic overview of CLASS4GL](https://class4gl.eu//wp-content/uploads/2019/01/image4-1024x794.png)

A CLASS4GL data package is available that can be directly used to perform and validate ABL model simulations and sensitivity experiments. The locations of the balloon soundings are performed for different climate regions as shown on the map.

[150 stations from IGRA of the reference dataset to perform and validate the ABL model simulations with CLASS4GL (see Sect. 2.2 of the CLASS4GL manuscript). The different climate classes are indicated with the colors according to the Köppen-Geiger climate classification. The markers indicate the locations of the atmospheric profiles from three observation campaigns (ie., HUMPPA, BLLAST and GOAMAZON)](https://class4gl.eu/wp-content/uploads/2019/01/image-1-480x300.png)]

[Data library of CLASS4GL](https://class4gl.eu/wp-content/uploads/2019/01/image-5-768x492.png)

### Reference
H. Wouters, I. Y. Petrova, C. C. van Heerwaarden, J. Vilà-Guerau de Arellano, A. J. Teuling, J. A. Santanello, V. Meulenberg, D. G. Miralles. A novel framework to investigate atmospheric boundary layer dynamics from balloon soundings worldwide: CLASS4GL v1.0. In preparation.


# Get started: 
see https://class4gl.eu/#getstarted


