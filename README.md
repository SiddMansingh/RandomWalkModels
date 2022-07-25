# Beginning of Website Making
## Exact Diagonalization Prerequisites
### Creation of Basis States and Lin Tables for easy reference
Before creating operations of the Hamiltonian, one needs to create basis states. The most crucial part is to find the index of a given configuration with the least time complexity, thus making the operation of a Hamiltonian on the basis state easier. A naive approach would be to evaluate $\sum_{i=0}^N 2^{\sigma_i}$ where $N$ is the number of spins at $\sigma_i= \pm 1$ representing the spin at site $i$. However one might need 
