## SOAP Compression

The length of the standard SOAP vector (powerspectrum)

$$p^{\alpha\beta}_{nn'l} = \sum_m c^\alpha_{nlm} c^\beta_{n'lm} = \bf{c}^\alpha_{nl} \cdot  \bf{c}^\beta_{n'l}$$

 is $\frac{1}{2}NS(NS+1)(L+1) + 1$ where $N$=`n_max`, $S$=`n_species`, $L$=`l_max` and the final element is `covariance_sigma0` (set after normalisation). Various compression strategies are available to reduce this $\mathcal{O}(N^2S^2)$ scaling.

## Tensor-reduction

The tensor-reduction introudced in [Tensor-reduced atomic density representations](https://doi.org/10.48550/arXiv.2210.01705) can be applied by first selecting which "density channels" (radial, species or both) are mixed together before selecting how the mixed channels should be coupled together (tensor product, element-wise). Tensor decomposition using radial-species mixing with element-wise coupling is written as

$$p_{kl} = \sum_m c_{klm} c_{klm} = \bf{c}_{kl} \cdot \bf{c}_{kl} $$

where
 $$c_{klm} = \sum_{\alpha n} W^k_{\alpha n} c^\alpha_{nlm} $$


| Keyword       | Values        | Description |
| -----------   | ------------- | ----------- |
| `R_mix`       | `T` or `F`    | mixes the radial channels |
| `Z_mix`       | `T` or `F`    | mixes the species channels |
| `sym_mix`     | `T` or `F`    | `T` means use "tensor-decomposition" whereas `F` means use "tensor-sketching" as described in the original [article](https://doi.org/10.48550/arXiv.2210.01705) |
| `K`           | int > 0       | How many mixed channels to create |
| `coupling`    | `T` or `F`    | `T` means use tensor product coupling across mixed channels, `F` means use element-wise|

Note that full tensor product coupling is always used across any "un-mixed" channels.

Examples
- `R_mix=T Z_mix=T K=5 sym_mix=T coupling=F` Means randomly mix the radial and species channels together to form 5 new channels then couple each of these channels to itself only. Length of final vector will be $K(L+1)$. If `coupling=T` then length would be $\frac{1}{2}K(K+1)(L+1)$ as each channel would be coupled to every other channel (per `l`).
- `R_mix=F Z_mix=T K=5 sym_mix=F coupling=F` Means randomly mix the species channels together so that there are `n_max*K` channels. Final length will be $\frac{1}{2}N(N+1)K(L+1)$ from element-wise coupling between the `K` mixed species channels and full tensor-product coupling between the radial channels.

## Radial and species sensitive correlation orders

The radial and species sensitive correlation orders can be set independently using

| Keyword       | Values | Description |
| -----------   | -------- | ----------- |
| `nu_R`       | `0, 1, 2`    | radially sensitive correlation order |
| `nu_S`       | `0, 1, 2`   | species senstitive correlation order |

The schematic below illustrates the physical interpretation and includes a correspondance between the `nu_R nu_S` notation used here and that used in the original article [Compressing local atomic neighbouhood descriptors ](https://www.nature.com/articles/s41524-022-00847-y).

<p align="center">
<img src="./translation_table.png" width="500">
</p>

## Other

- Setting `diagonal_radial=T` only includes terms where $n=n'$ so that the length becomes $\frac{1}{2}NS(S+1)(L+1) + 1$
- Different species can be grouped together using the `Z_map` keyword where commas separate groups of elements e.g. `Z_map={8, 23 41 42}` treats all the metals as identical and Oxygen as distinct. As the power spectrum is correlation order 2 it is also possible to specify two distinct groupings separated by a colon e.g. `Z_map={8, 23, 41, 42: 8, 23 41 42}` where in the first density every element is distinct whilst in the second again the metals are all treated as identical but Oxygen is distinct.
