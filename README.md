# Beginning of Website Making
## Exact Diagonalization Prerequisites
### Creation of Basis States and Lin Tables for easy reference
Before creating operations of the Hamiltonian, one needs to create basis states. The most crucial part is to find the index of a given configuration with the least time complexity, thus making the operation of a Hamiltonian on the basis state easier. A naive approach to find this mapping between a state and it's index would be to evaluate $\sum^N_{i=0} \sigma_i\cdot 2^{i}$ where $N$ is the number of spins and $\sigma_i= \pm 1$ represents the spin at site $i$. This expression returns a unique index for each configuration, however the evaluation of  $2^i$ might be limited by memory. 

[Lin lookup tables](https://aip.scitation.org/doi/pdf/10.1063/1.4823192) offer an easier solution to this hashing problem. The inclusion of symmetries of the Hamiltonian reduces the size of the basis and hence the memory occupied by the tables. The code demonstrates the generation of basis states for a $n$ spin system with $S_z$ symmetry. One can feed in the arguments $n=$ Size of the system and $k=$ Number of up spins for which the basis would be created.


In the upcoming weekend, the implementation for Exact Diagonalization would be made public.
