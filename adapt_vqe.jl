using MPI
MPI.Init()
const root = 0
const tag_request = 1
const tag_idx = 2
using LinearAlgebra
push!(LOAD_PATH, "./QuantumSpins/src")
push!(LOAD_PATH, "./QuantumCircuits/src")
push!(LOAD_PATH, "./MPSSimulator/src")
push!(LOAD_PATH, "./DMV/src")
using QuantumSpins
using QuantumCircuits
using MPSSimulator
using DMV
using Zygote
using NLopt

if length(ARGS)>1
    D_global = parse(Int, ARGS[2])
else
    D_global = 128
end
if length(ARGS)>2
    Group_size = parse(Int, ARGS[3])
else
    Group_size = 4
end
if MPI.Comm_rank(MPI.COMM_WORLD) == root
	if D_global <= 0
		error("Bond dimension must > 0!")
	end
	if Group_size <= 0
		error("Group size must > 0!")
	end
    println("Bond dimension = $D_global")
	println("Group size = $Group_size")
end





function vqe(n_qubits::Int, n_params::Int, map::Vector{Tuple{Float64, Int}}, circuit::QCircuit, observables::Vector{Hamiltonian}, comm :: MPI.Comm ;
                init_params::Vector{Float64}=[], initial_state::MPS, trunc::MPSTruncation)

	rank = MPI.Comm_rank(comm)
	size = MPI.Comm_size(comm)
	(init_params == []) && (init_params = zeros(Float64, n_params))
	num_groups = length(observables)
	

    if rank == root
	        
            
            function fmw(x::Vector, grad::Vector)
				MPI.bcast(x, root, comm)
				  #println("At $(x)")
			
						n_task = length(observables)
						idx = 1
						recv = zeros(1)

						req = MPI.Irecv!(recv, MPI.ANY_SOURCE, tag_request, comm)
						while (n_task > 0)
							flag, status = MPI.Test!(req)
							if flag
								MPI.isend(idx, status.source, tag_idx, comm)
								idx += 1
								n_task -= 1
								req = MPI.Irecv!(recv, MPI.ANY_SOURCE, tag_request, comm)
							end
						end
					
						n_worker = size - 1
						while (n_worker > 0)
							flag, status = MPI.Test!(req)
							if flag
								MPI.isend(-1, status.source, tag_idx, comm)
								n_worker -= 1
								if n_worker > 0
									req = MPI.Irecv!(recv, MPI.ANY_SOURCE, tag_request, comm)
								end
							end
						end
			

			buffer = MPI.Reduce(zeros(n_params+1), +, root, comm)

					v = buffer[1]
					grad[:] = buffer[2:end]
				  println("got $v, ||g|| = $(norm(grad))"); flush(stdout)
				v
			end

			opt = Opt(:LD_LBFGS, n_params)
            opt.min_objective = fmw
            opt.maxeval = 50
            opt.ftol_abs = 1e-6
            opt.vector_storage = opt.maxeval + 2

            (minf, minx, ret) = optimize(opt, init_params)
			MPI.bcast(nothing, root, comm)
            MPI.bcast((minf, minx), root, comm)
    else
		while (true)
			x0::Union{Nothing, Vector{Float64}} = MPI.bcast(nothing, root, comm)
			if x0 == nothing
				break
			end
			datapack_local = zeros(Float64, n_params + 1)

			while (true)
				MPI.Isend(0, root, tag_request, comm)
				while ((flag, idx, status) = MPI.irecv(root, tag_idx, comm))[1] == false end

				if idx == -1
					MPI.Reduce(datapack_local, +, root, comm)
					break
				end


				_datapack, _t = get_expectation_and_gradient(x0, map, circuit, initial_state, observables[idx], trunc, timing=false)
				datapack_local += _datapack

			end
		end	
		minf, minx = MPI.bcast(nothing, root, comm)
    end

	circuit2 = copy(circuit)
	QuantumCircuits.reset_parameters!(circuit2, DMV.apply_map(minx, map=map))
	psi = MPSSimulator.apply(circuit, initial_state, trunc=trunc)
	QuantumSpins.normalize!(psi)
	return minf, minx, psi
end




function adapt_vqe(parameters :: Dict, comm :: MPI.Comm; adapt_max_iter = 50,
                                                         adapt_g_tol = 1e-3,
                                                         adapt_f_tol = 1e-7,
							 batch_size = 1,
                                                         trunc = MPSTruncation(ϵ=1.0e-6, D=D_global))

	rank = MPI.Comm_rank(comm)
	n_procs = MPI.Comm_size(comm)
	n_qubits::Int = parameters["n_qubits"]
	n_operators_pool::Int = parameters["n_params"]
	spin_orbital_occupied_indices::Vector{Int} = parameters["spin_orbital_occupied_indices"]
	operator_pool::Vector{QubitOperator} = parameters["ucc_operator_pool_qubit_op"]	
	observables::Vector{Hamiltonian} = parameters["hamiltonian_qubit_op"]
	whole_ham::Hamiltonian = parameters["whole_ham"]
	alg = StableArith(D=trunc.D, ϵ=trunc.ϵ)
	psi = hartreefock_state_mps(n_qubits, spin_orbital_occupied_indices, trunc)
	hf_state = copy(psi)

    if rank == root
        println("qubits: ", n_qubits)
        println("parameters: ",length(operator_pool))
        println("global hamiltonian terms: ",whole_ham.data.count)
	println("Building MPOs of operator pool")
    end
	

	# Build MPOs of operator pool
	operator_pool_temp = Hamiltonian.(operator_pool)
	operator_pool_mpos = MPSSimulator._MPO.(length(hf_state), operator_pool_temp)

	operators = Vector{QubitOperator}()
	operator_indices = Vector{Int}()
	θs = Vector{Float64}()
    opt_energy, opt_energy_last = 0.0, 0.0
    
t1=time_ns()
    for iter in 1:adapt_max_iter
       
        if rank == root
            println("=========ADAPT-VQE step $(iter)=============")
        end


		# ================= Measure residual gradients using expectation=========
		gs_local = zeros(Float64, n_operators_pool)
		for i in p_range(n_operators_pool, n_procs, rank + 1)
			Op_psi, err = mpompsmult(operator_pool_mpos[i], psi, alg)
			gs_local[i] = abs(2*real(expectation(psi, whole_ham, Op_psi)))
		end
		gs = MPI.Allreduce(gs_local, +, comm)
		# =======================================================================


        # ================ Check norm of residual gradient ======================
		g_norm = norm(gs)
        if rank == root
            println("Residual gradient norm = $(g_norm)")
        end

        if g_norm <= adapt_g_tol
            if rank == root
                println("Norm of residual gradient = $(g_norm) <= $(adapt_g_tol), ADPAT-VQE completed.")
            end
            break
        end
        # =======================================================================


        # ================ Batched add operators ================================
		op_indices = sortperm(gs, rev=true)[1:batch_size]
		for idx in op_indices
            push!(operators, operator_pool[idx])
			push!(operator_indices, idx)
			push!(θs, zero(Float64))
        end
        if rank == root
            println("Add $(op_indices) -th operators to the ansatz")
			println("Operator pool: $(operator_indices)")
        end
        # =======================================================================



        # ================ Grow ansatz and perform VQE ==========================
        vqe_circuit, vqe_map = construct_vqe_circuit(operators)
        opt_energy, θs, psi = vqe(n_qubits, length(θs), vqe_map, vqe_circuit, observables, comm, 
				init_params=θs, initial_state=hf_state, trunc=trunc)

		# =======================================================================
        

        # ================ Check the adapt_f_tol ================================
        if rank == root
            println("Optimized energy: $(opt_energy)")
        end
        if (iter > 1) && (abs(opt_energy - opt_energy_last) <= adapt_f_tol)
            if rank == root
                println("The energy difference between two iterations = |$(opt_energy) - $(opt_energy_last)| <= $(adapt_f_tol), \
                ADPAT-VQE completed.")
            end
            break
        end
        opt_energy_last = opt_energy
        # =======================================================================
         
        if rank == root
            println("==========================================")
        end
        
    end
t2=time_ns()
    if rank == root
        println("ADPAT-VQE result: $(opt_energy), total time: $((t2-t1)/1e9)")
    end
end


function main()

	comm = MPI.COMM_WORLD
	rank = MPI.Comm_rank(comm)
	p = Dict{String, Any}()
	if rank == root
		f = open(ARGS[1], "r")
		p = read_binary_dict(f)


		q = read_binary_qubit_op(f, p["n_qubits"])
		println("original global hamiltonian terms: ",length(q))


		H = Hamiltonian(q)
		println("global hamiltonian terms: ",H.data.count)
        p["whole_ham"] = H



			p["hamiltonian_qubit_op"] = make_groups(H, p_range(length(H.data), 1, 1), group_size=Group_size)
			MPI.bcast(p, root, comm)

	else
			p = MPI.bcast(nothing, root, comm)
	end

	adapt_vqe(p, comm)
end
main()
