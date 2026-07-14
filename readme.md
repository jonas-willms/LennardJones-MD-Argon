## General info

An example project to plot the aggregate states of "argon" and their neighbor distributions.

It is inspired by the orignal paper:

"Computer "Exyeriments" on Classical Fluids. I. Thermodynamical Properties of Lennard-Jones Molecules"
by Loup Verlet 1967

The example script simulates "argon" atoms at 3 different density and temperature combinations.
The result will be a figure showing the three states ().

## How to Run

To run the simulation and generate the phase distribution plot, you need to execute the script using the local project environment.

Open your terminal, navigate to the root directory of this repository (`ArgonMD/`), and run the following command:

```bash
julia --project=. scripts/plot_phases.jl
```
