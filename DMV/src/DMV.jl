module DMV

using MPI
using LinearAlgebra
using QuantumSpins
using QuantumCircuits
using MPSSimulator
using Zygote
using Random

export MPSTruncation


export QubitOperator, Hamiltonian, flattening!, bond_dimension, make_groups
include("hamiltonian.jl")


export p_range, hartreefock_state_mps, xzqsize, stdev, mean
include("utils.jl")

export read_binary_qubit_op, read_binary_dict
include("io.jl")

include("communication.jl")

export construct_vqe_circuit, construct_hea_circuit
include("circuit.jl")


export get_expectation, get_expectation_and_gradient
include("vqe.jl")

end