function read_binary_qubit_op(f::IO, n_qubits::Int) :: QubitOperator

    pauli_symbol_dict = Dict(
        0 => "I",
        1 => "X",
        2 => "Y",
        3 => "Z"
    )
	qubit_op_dict = QubitOperator()
	len = read(f, Int64)
	if len > 0			# paulis are stored in dense array
		for i in 1:len
			coeff_tmp = read(f, ComplexF64)
			pauli_str_tmp = Int8[read(f, Int8) for i in 1:n_qubits]
			pauli_str_tuple = Tuple([(i, pauli_symbol_dict[pauli_str_tmp[i]])
                                 for i in 1:n_qubits
                                 if pauli_str_tmp[i] != 0])
			push!(qubit_op_dict, pauli_str_tuple => coeff_tmp)
		end
	else				# paulis are stored in compressed array
		len = -len
		for i in 1:len
			coeff_tmp = read(f, ComplexF64)
			len_pauli_array = read(f, Int32)
			pauli_str_tuple = Vector{Tuple{Int64, String}}()
			for i in 1:len_pauli_array
                # Indices are shifted by 1 due to array indices start at 1 in Julia.
				push!(pauli_str_tuple, (read(f, Int32) + 1, pauli_symbol_dict[read(f, Int8)]))
			end
			push!(qubit_op_dict, Tuple(pauli_str_tuple) => coeff_tmp)
		end
	end
	return qubit_op_dict
end

function read_binary_dict(f::IO)::Dict{String,Any}

    identifier = read(f, Float64)
    err_msg = identifier != Float64(99.0212) && error("The file is not saved as vqe parameters.")

	d = Dict{String, Any}()
	d["n_qubits"] = Int64(read(f, Int32))
	len_spin_indices = read(f, Int32)

    # Indices are shifted by 1 due to array indices start at 1 in Julia.
	d["spin_orbital_occupied_indices"] = Int64[read(f,Int32) + 1 for i in 1:len_spin_indices]

	d["n_params"] = Int64(read(f, Int32))
	d["ucc_operator_pool_qubit_op"] = Vector{QubitOperator}(undef, d["n_params"])
	for i in 1:d["n_params"]
		d["ucc_operator_pool_qubit_op"][i] = read_binary_qubit_op(f, d["n_qubits"])
	end
	return d
end