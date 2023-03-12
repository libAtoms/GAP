.. _accelerated-aimd:

*******************************************************
Acceleration of Ab-initio MD with On-the-Fly GAP Models
*******************************************************

This is the documentation of the ``hybrid-md`` package bundled with GAP, which   provides a general framework of 
accelerating Ab-initio MD (AIMD) with Machine Learning models generated "On-the-Fly".

The framework is intended to be agnostic of the AIMD code and the ML code as well. For the time being, an interface 
between the CASTEP code and GAP is available publicly, integration of other ab-initio codes and ML frameworks is in 
the making. If you would have any software/research ideas, or would like to contribute in any way, please feel free 
to reach out to with G. Csányi and T. K. Stenczel for further details.

Getting started & Usage
***********************

.. toctree::
    :maxdepth: 2

    accelerated-aimd_installation.rst
    accelerated-aimd-input-reference.rst
    accelerated-aimd_decision-making.rst


Program structures and design
*****************************

The program is intended to be user-accessible and extendable. CASTEP implements
an interface to external force evaluation codes - including QUIP/GAP.

Design Principles
#################

The interface provides acceleration of an AIMD simulation, according to the following design principles:

1. MD is driven by the AIMD code (this is the main program)

2. Decision making program independent of AIMD code

3. ML code is used by the AIMD code only for force evaluation, and by decision making code for updating of models

4. The system is agnostic of programming languages used, the decision making code could be implemented in any other language as long as the API is the same

Extendability & Development
***************************

Users can extend the framework in multiple ways:

#. Build custom decision making routines: e.g.
    - system specific adaptations
    - utilising confidence estimates of models
    - any external tools working out when the model is extrapolating / could use more data: e.g. sample more when connectivity changes - bond break/form
#. Support other AIMD codes: There are no specific requirements apart from an internal MD loop being available. Please get in touch with G. Csányi and T. K. Stenczel for further details.
#. Support for other ML tools - the framework is not depending on GAP

