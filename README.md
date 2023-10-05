# Large-scale heterogeneities can alter the characteristics of compressive failure and accelerated seismic release
[DOI: https://doi.org/10.1103/PhysRevE.108.014131](https://doi.org/10.1103/PhysRevE.108.014131)

Authors: 

Drew Patton      -- Code, Analysis, and writer

Joern Davidsen   --      Supervisor and writer

Thomas Goebel    --    Experimenter and writer

Grzegorz Kwiatek --    Experimenter and writer

## Abstract
Externally stressed brittle rocks fail once the stress is sufficiently high. This failure is typically preceded by
a pronounced increase in the total energy of acoustic emission (AE) events, the so-called accelerated seismic
release. Yet, other characteristics of approaching the failure point such as the presence or absence of variations
in the AE size distribution and, similarly, whether the failure point can be interpreted as a critical point in a
statistical physics sense differs across experiments. Here, we show that large-scale stress heterogeneities induced
by a notch fundamentally change the characteristics of the failure point in triaxial compression experiments under
a constant displacement rate on Westerly granite samples. Specifically, we observe accelerated seismic release
without a critical point and no change in power-law exponent $\epsilon$ of the AE size distribution. This is in contrast to
intact samples, which exhibit a significant decrease in $\epsilon$ before failure. Our findings imply that the presence or
absence of large-scale heterogeneities play a significant role in our ability to predict compressive failure in rock.

## Project Description
My research delves into the intricate behavior of rocks when subjected to an increase in significant stress, such as that experienced during geological events like earthquakes.
The process by which rocks fracture under stress is far from straightforward; instead, it is quite a complex process. Imagine a piece of rock with tiny cracks within it. 
As external pressure mounts, these cracks form, grow, and coalesce, releasing stored energy in a manner analogous to the distinct crackle noise when stepping on dry leaves.

What I'm undertaking is an investigation into how the presence of a pre-existing notch in the rock's structure influences the way it breaks apart.
To accomplish this, we employed specialized equipment and instrumentation to record the acoustic emissions produced by the rock as it undergoes stress. 
These acoustic signals provide valuable insights into the internal processes at play. In particular, acoustic emissions provide a wealth of knowledge in the number of microcrack events over a period of time, the acoustic energy released, the position in 3-dimensions,
and even discern how these microcracks form – whether through compaction, shearing, or tensile forces.
The rock samples we used are all Westerly granite samples, and they are separated into two categories: notched and intact.

My findings reveal that the fracture patterns in rocks can vary substantially depending on their structural properties: Intact samples follow a failure process known as a critical failure, while the notched samples do not. This understanding carries significant implications for the prediction and mitigation of geological hazards like earthquakes and landslides, ultimately contributing to enhanced safety measures.

## Folder Description

### RawData 
Contains the raw data obtained/received by the triaxial compression experiments.
Each folder in raw data is a Westerly granite sample used in the analysis. Samples with an 'n' are the notched samples, the remaining are intact samples.
Each sample folder contains two different datasets labeled: 'AE' for the acoustic emission data and 'MTS' for the force-time data.
The dataset folders may either contain a .mat file or .txt file (this was how they were received).

For the datasets 'AE,' the only data of interest is the time: 'Time', the adjusted amplitde: 'AdjAmp' for the energy, the 3-dimensional coordinates: 'x', 'y', and 'z', and the polarity: 'pol'.

For the datasets 'MTS', the only data of interest is the time: 'Time', and the force: 'Force' or 'force'.

### External_Code
Contains fortran code in the form of the file: '2D_EM_continuum_p_val.f90' for the quick application of a 'Maximum likelihood estimate' method to determine the best power law exponent fit.

The 'a.exe' file executes the fortran code from Jupyter.

The only file of interest is the '2D_EM_continuum_p_val.f90' and the output file that ends with extension: '.EMap'. The remaining may be ignored for this project.

## Main File


