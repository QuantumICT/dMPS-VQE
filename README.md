# dMPS-VQE

## Requirement
### Julia software
* Julia 1.7+

### Julia Packages
* MPI v0.20.5+ (with mpiexecjl installed)
* NLopt v0.6.5+
* Parameters v0.12.3+
* TensorOperations v3.2.4+
* Zygote v0.6.51+
* KrylovKit v0.6.0+

## How to run
For example, to run ADAPT-VQE on H$_2$ molecule with 4 MPI processes, simply run
```
mpiexecjl -n 4 julia adapt_vqe.jl inputs/input_h2
```