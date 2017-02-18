[![Build Status](https://travis-ci.org/HomerReid/GLMTSolver.svg?branch=master)](https://travis-ci.org/HomerReid/GLMTSolver)

# GLMTSolver
========

This is a simple code that implements a generalized Lorenz-Mie theory (GLMT)
solver for electromagnetic scattering from spherically-symmetric bodies.

The geometries that may be handled consist of arbitrarily many
spherical layers, of arbitrary inner and outer radii, each with
arbitrary frequency-dependent permittivity.

The allowed incident fields are

 + plane waves

 + spherical waves

The outputs that may be computed include 

 + scattered and total **E** and **H** fields at arbitrary points in space

 + scattered and absorbed power, force and torque as obtained by integrating the Poynting vector and Maxwell stress tensor over a bounding sphere surrounding the particle

# Theory
========

The theoretical approach and notation used by this code 
is described in [this memo][scuffSpherical]

# Specifying geometries
========

Geometries are described by simple text files conventionally given
file extension `.GLMT.` Blank lines and comments (lines starting with `#`)
are skipped. Each line specifies a single spherical layer, with
the following syntax

```
 RMAX   MATERIAL
```

where 

 + `RMAX` is the outer radius of the layer in microns

 + `MATERIAL` is a [<span font-variant="small-caps">scuff-em</span> material designation][scuffMaterials].

Here's a `.GLMT` file describing a sphere of radius 1 &mu;m with dielectric constant 10:

````
1.0	CONST_EPS_10
````

And here's a `.GLMT` file describing the
 [doubly-resonant multilayered spherical particle](WadePaper)
considered by C. W. Hsu et al, "Theoretical criteria for
scattering dark states in nanostructured particles:"

````
0.100 	SiO2
0.051 	Silver
0.031 	SiO2
0.020	Silver

````

If one of the layered material regions in your geometry is
actually a void (no material), give it material property VACUUM.
For example, here's a geometry describing the layered particle
above, but with the silver removed to yield just two concentric 
spherical shells of SiO2.

````
0.100 	SiO2
0.051 	VACUUM
0.031 	SiO2
0.020	VACUUM

````

# Structure of the code
========

`GLMT-scatter.cc` builds to yield a standalone executable
command-line code that implements a subset of the functionality
of [<span font-variant="small-caps">scuff-scatter</span>][scuffScatter]
(see examples below).

The main solver is implemented in the code file `GLMTSolver.cc`
The routine `CreateGLMTSolver()` parses a `.GLMT` file and constructs
a data structure named `GLMTSolver` containing all information
the scattering geometry. This structure contains tables of 
spherical-wave coefficients which are initialized for a 
given frequency and a given incident field by the routine
`Solve()'. Once `Solve()` has been called, you can 
call post-processing routines:

 + `GetFields()` computes the scattered and total **E** and **H** fields
   at arbitrary points in space

 + `GetDSIPFT()` computes the power, force, and torque on the body,
using the **D**isplaced **S**urface **I**ntegral computational strategy
described [in this paper][PFTPaper].

# Examples
========

## Power scattered by doubly-resonant
L=2

ARGS=""
ARGS="${ARGS} --GLMTFile ${GLMTDIR}/WadeParticle4.GLMT"
ARGS="${ARGS} --PWDirection 0 0 1"
ARGS="${ARGS} --PWPolarization 1 0 0"
ARGS="${ARGS} --Omega 12.5"
ARGS="${ARGS} --EPFile ${EPDIR}/TP45"
ARGS="${ARGS} --FileBase WadeParticle4.GLMT.L${L}"
ARGS="${ARGS} --LMax ${L}"
${PREFIX} ${DIR}/${CODE} ${ARGS}


##################################################
for L in 1 2 3 
do
ARGS=""
ARGS="${ARGS} --GLMTFile ${GLMTDIR}/WadeParticle4.GLMT"
ARGS="${ARGS} --PWDirection 0 0 1"
ARGS="${ARGS} --PWPolarization 1 0 0"
ARGS="${ARGS} --OmegaFile L100800"
ARGS="${ARGS} --PFTFile WadeParticle4.L${L}.GLMTPFT"
ARGS="${ARGS} --LMax ${L}"
ARGS="${ARGS} --DSIRadius 0.2"
${PREFIX} ${DIR}/${CODE} ${ARGS}

scuffMaterials:		http://homerreid.github.io/scuff-em-documentation/reference/Materials
scuffSpherical:		http://homerreid.github.io/scuff-em-documentation/tex/scuffSpherical.pdf
