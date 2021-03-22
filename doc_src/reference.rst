**********************
SOAP and GAP reference
**********************

This page collects useful information in working with the GAP code. 

SOAP vectors
************

The description of the local environment by SOAP vectors is as follows. The neighbour density around each atom is expanded in a set of radial and angular basis functions, with a different expansion for each atomic species in the neighbour environment. The SOAP vector is then formed from these expansion coefficients. The array has two radial and one angular indices. The number of radial basis functions is `n_max`, the number of angular basis functions `l_max+1`. The order of the elements in the vector is built as follows.

```
for i=1:n_species*n_max
  for j=i:n_species*n_max
    for l=0:l_max
```


