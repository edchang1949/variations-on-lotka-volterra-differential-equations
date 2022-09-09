# variations-on-lotka-volterra-differential-equations

Some code aiming to simulate the Lotka-Volterra predator-prey differential equations and to resolve some assumptions and limitations in the original equations using Python and Matplotlib.

The original equations are from https://web.ma.utexas.edu/users/davis/375/popecol/lec10/lotka.html

## Overview
There are three parts to the code. It is recommended that only one part is run at any given time, so the second and third parts are commented out to start.

1. A simulation of the original equations using 4th order Runge-Kutta approximation.

2. This part features new parameters that add needed complexity to the original equations. There are now two predator species and two prey species. The two predator species can interact with each other such that overcrowding causes their populations to decrease, and similarly, too many prey overall will result in decreasing prey populations. The strength of these interactions can be tuned using the interaction coefficients "hic" for prey and "pic" for predators. Setting these to zero will essentially cause the predator and prey species to ignore the population of their peers. Additionally two more parameters have been added to allow each predator species to prefer one prey over another. This can be tuned using the cross coefficients "cc12" between prey one and predator two and "cc21" between prey two and predator one. Setting these to zero will have predator one only affecting prey one and predator two only affecting prey two.

3. This part aims to fix a shortcoming of part two by allowing the population of each species to be capped at certain population sizes. This adds more realism to the somewhat wild simulations that can come from part two.