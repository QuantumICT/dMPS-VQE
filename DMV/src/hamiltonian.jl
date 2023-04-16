const QubitOperator = Vector{Pair{Tuple{Vararg{Tuple{Int64, String}}}, ComplexF64}}
const Hamiltonian = QuantumCircuits.QubitsOperator
function Hamiltonian(x::QubitOperator)
	qop = Hamiltonian()
	#if length(x) == 0
	#	QuantumCircuits.add!(qop, QubitsTerm(Dict(1=>"I"), coeff=0.0))
	#end
	for (k, v) in x
		d = Dict{Int, Union{AbstractString, AbstractMatrix}}()
		for (idx, qubit_gate) in k
			d[idx] = qubit_gate
		end
		if d.count == 0
			d[1] = "I"
		end
		QuantumCircuits.add!(qop, QubitsTerm(d, coeff=v))
	end
	return qop	
end

function make_groups(QO::Hamiltonian, r::UnitRange{Int}; group_size::Int)::Vector{Hamiltonian}
	g = Vector{UnitRange{Int}}()
	pairs = collect(QO.data)
	for i in r.start : group_size : r.stop
		push!(g, i:min(i+group_size-1,r.stop))
	end
	qo = [Hamiltonian(QuantumCircuits.QOP_DATA_TYPE(@view pairs[r])) for r in g]
	if length(qo) == 0
		qo = [Hamiltonian(QubitOperator())]
	end
	return qo
end

function randomly_swap_terms!(H1::Hamiltonian, H2::Hamiltonian)
	k1 = rand(keys(H1))
	k2 = rand(keys(H2))
	v1 = copy(H1.data[k1])
	v2 = copy(H2.data[k2])
	delete!(H1.data, k1)
	delete!(H2.data, k2)
	nv1 = get!(H1.data, k2, QuantumCircuits.QOP_DATA_VALUE_TYPE())
	nv2 = get!(H2.data, k1, QuantumCircuits.QOP_DATA_VALUE_TYPE())
	append!(nv1, v2)
	append!(nv2, v1)
end

function bond_dimension(n_qubits::Int, H::Hamiltonian)
    return QuantumSpins.bond_dimension(MPSSimulator._MPO(n_qubits, H))
end


function flattening!(n_qubits::Int, Hs::Vector{Hamiltonian}, bds::Vector)
    if length(bds) < 1
        resize!(bds, length(Hs))
        bds[:] = bond_dimension.(n_qubits, Hs)
    end
    (length(bds) <= 1) && (return)
    imin, imax = argmin(bds), argmax(bds)
    (imin == imax) && (return)
    randomly_swap_terms!(Hs[imin], Hs[imax])
    bds[imin] = bond_dimension(n_qubits, Hs[imin])
    bds[imax] = bond_dimension(n_qubits, Hs[imax])
end