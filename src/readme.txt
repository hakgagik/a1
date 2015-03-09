Gagik Hakobyan
gh328
CS 5143 Assignment 1

Implementation details:

My implementation closely followed the procedure given in lecture and in the Macklin and Muller paper[1]. I've discussed the project with Jack and Adam (Jack is the one that reminded me to make sure particles don't interact through thin walls in my artifact). 


The acceleration structure I'm using is the one suggested. The scene is split into cubic cells with side length h. At the start of the simulation loop, each particle is assigned a cell. When a particle needs to check its neighbors, I check its own cell, as well as its 26-connected neighboring cells for particles closer than h.

The main differences between my implementation and the Macklin paper are in vorticity confinement and viscosity. When we blur anything using the blurring kernels (poly6 and spiky), we generally multiply by m/rho, but the Macklin paper doesn't seem to do that when calculating vorticity or the viscosity correction. If this is omitted, the sums themselves yield VERY large values, meaning c and epsilon must be very small (much smaller than c = 0.01 as suggested). So I made a point to divide by the local density when using either of the blurring kernels. As a result, my vorticity implementation is closer to the one given by Muler, et al [2].

I used c = 0.05 and epsilon_viscosity = 2. This seemed to give the optimum combination of splashing and cohesiveness. All of the other parameters were kept the same as those suggested on piazza.

Performance:

My implementation runs smoothly (30+ fps) when under ~500 particles, such as with cube_drop.txt. It runs acceptably (~10 fps) at around 2500 particles. I believe the implementation could be made faster by storing particle distances with their neighbors (since distance calculations are the second most common operations, second only to kernel evaluations).

Artifact: Leaky washing machine

The artifact is inside a 4x2x2 simulation cell. Unfortunately, I don't remember how many particles I'm using, but it should be on the order of 12000. There is a thin wall at x = 2, with a slit in the wall between z=0 and z=0.2. Particles that collide with the wall have their velocity in the normal direction set to 0 (by having a force proportional to their velocity in the normal direction applied to them) and their deltaP is changed to land the particle an epsilon away from the wall. (I account for the wall in particle interactions by making sure the segment between the two particles doesn't intersect the wall.)

The artifact (artifact.webm) starts with dam break, then has more particles added to it over time in clumps. The fluid first splashes into the wall, then slowly leaks out through the slit into the other half of the simulation cell. The washing machine (a rotational force field with f=(-rsin(theta), 0, rcos(theta)), where r is distance from the center of the cell and theta is the angle from the z-axis) turns on at 3 seconds, spinning the fluid around and forcing it through the slit. The artifact demonstrates splashing very well. I especially like how the fluid, spreading after it flows out of the slit, causes a swirling effect near the slit. 

The artifact was rendered with the Masuda renderer, then compiled into a webm file at 60 fps. 

I also included an artifact I created yesterday in a 1x1x1 cell, but was not very satisfied with (my.webm).


[1] Miles Macklin and Matthias Müller. Position Based Fluids. ACM Trans. Graph. 32, 4, Article 104 (July 2013), 12 pages

[2] Matthias Müller, David Charypar, Markus Gross, Particle-based fluid simulation for interactive applications, 2003 ACM SIGGRAPH / Eurographics Symposium on Computer Animation (SCA 2003), August 2003, pp. 154-159