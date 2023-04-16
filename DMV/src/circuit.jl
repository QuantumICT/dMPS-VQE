function construct_hea_circuit(n_qubits::Int, n_layers::Int)
    circuit = QCircuit()
    m = Vector{Tuple{Float64, Int}}()
    param_idx = 0
    for n_l in 1:n_layers
        for qi in 1:n_qubits
            push!(circuit, RyGate(qi, 0.0, isparas=true))
            param_idx += 1
            push!(m, (1.0, param_idx))
        end
        for qi in 1:n_qubits
            push!(circuit, RzGate(qi, 0.0, isparas=true))
            param_idx += 1
            push!(m, (1.0, param_idx))
        end

        # A circular entanglement as:
        # CNOT(1, 2), CNOT(2, 3), ..., CNOT(N-1, N), CNOT(N, 1)
        for qi in 1:n_qubits-1
            push!(circuit, CNOTGate(qi, qi+1))
        end
        push!(circuit, CNOTGate(n_qubits, 1))
    end

    # Add a last layer of single qubit gates.
    for qi in 1:n_qubits
        push!(circuit, RyGate(qi, 0.0, isparas=true))
        param_idx += 1
        push!(m, (1.0, param_idx))
    end
    for qi in 1:n_qubits
        push!(circuit, RzGate(qi, 0.0, isparas=true))
        param_idx += 1
        push!(m, (1.0, param_idx))
    end

    return circuit, m
end

# Return circuit and map
function construct_vqe_circuit(operator_pool::Vector)
	circuit = QCircuit()
	m = Vector{Tuple{Float64, Int}}()
	for i in 1:length(operator_pool)

		for (terms,coeff) in operator_pool[i]
			length(terms) == 0 && continue
			for (idx, qubit_gate) in terms
				if qubit_gate == "X"
					push!(circuit, HGate(idx))
				elseif qubit_gate == "Y"
					push!(circuit, gate(idx, hy_matrix))
				end
			end

			for i in length(terms)-1:-1:1
				push!(circuit, CNOTGate(terms[i+1][1], terms[i][1]))
			end

			push!(circuit, RzGate(terms[1][1], 0.0, isparas=true))
			push!(m, (-2.0*imag(coeff),i))

			for i in 1:length(terms)-1
				push!(circuit, CNOTGate(terms[i+1][1], terms[i][1]))
			end

			for (idx, qubit_gate) in terms
				if qubit_gate == "X"
					push!(circuit, HGate(idx))
				elseif qubit_gate == "Y"
					push!(circuit, gate(idx, hy_matrix))
				end
			end

		end

	end
	return circuit, m
end

# Zygote does not treat keywords as variables
# Parameters after ";" are treated as keywords in Julia
function reset_parameters_wrapper(p::Vector{<:Real}; x::QCircuit) :: QCircuit
	reset_parameters!(x, p)
	return x
end
Zygote.@adjoint reset_parameters_wrapper(p::Vector{<:Real}; x::QCircuit) = reset_parameters_wrapper(p, x = x), z -> (z,)


# Apply the map to variational parameters in VQE
# The adjoint function is not necessary since derivative can be generated automatically by Zygote.jl
function apply_map(amps :: Vector{Float64}; map :: Vector{Tuple{Float64, Int}}) :: Vector{Float64}
	return [line[1]*amps[line[2]] for line in map]
end
