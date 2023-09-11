# The stabilizer group

We use a 2D Bool array to store generators of the stabilizer group.
The first index labels the stablizer; The second index expresses the element by 1-site Pauli operators.
    Note that for n^th d.o.f. (starting from 1), 2*n-1 is for Pauli_X, 2*n is for Pauli_Z.