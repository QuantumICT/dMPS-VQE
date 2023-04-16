const hy_matrix = [+1.0+0.0im +0.0-1.0im; 0.0+1.0im -1.0+0.0im] * 0.5 * 2^0.5

function p_range(n::Int, n_procs::Int, pid::Int)::UnitRange{Int}
	aprocs = n_procs - n % n_procs + 1
	q = n ÷ n_procs
	if (pid < aprocs)
			pstart = (pid - 1) * q + 1
			pend = pstart + q - 1
	else
			pstart = (aprocs-1) * q + (pid - aprocs)*(q+1) + 1
			pend = pstart + q
	end
	return pstart:pend
end

# Normalization of MPS
normalize_mps(psi::MPS) = normalize(psi)
Zygote.@adjoint normalize_mps(psi::MPS) = normalize_mps(psi), z -> (z,)

function hartreefock_state_mps(n_qubits::Int, spin_orbital_occupied_indices::Vector{Int}, trunc::MPSTruncation)::MPS
	hf_state = statevector_mps(n_qubits)
	apply!(QCircuit([XGate(i) for i in spin_orbital_occupied_indices]), hf_state, trunc=trunc)
	return hf_state
end

xzqsize(x) = Base.format_bytes(Base.summarysize(x))

stdev(a::Vector) = begin
	L = length(a)
	m = sum(a) / L
	sqrt(sum((a .- m).^2)/L)
end

mean(a::Vector) = begin
	L = length(a)
	sum(a) / L
end