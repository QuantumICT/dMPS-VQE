function MPI.Scatterv(H :: Union{Nothing, Hamiltonian}, root :: Integer, comm::MPI.Comm; group_size::Int) :: Vector{Hamiltonian}
	if root == MPI.Comm_rank(comm)
			@assert H != nothing
			len = H.data.count
			n_procs = MPI.Comm_size(comm)
			reqs = Vector{MPI.Request}(undef, n_procs - 1)
			for i in 1:n_procs-1
					hs = make_groups(H, p_range(len, n_procs, i+1), group_size=group_size)
					reqs[i] = MPI.isend(hs, i, 0, comm)
			end
			hs_self = make_groups(H, p_range(len, n_procs, 1), group_size=group_size)
			MPI.Waitall!(reqs)
			return hs_self
	else
			hs, status = MPI.recv(root, 0, comm)
			return hs
	end
end

function MPI.Scatterv(Hs :: Union{Vector{Nothing}, Vector{Hamiltonian}}, root :: Integer, comm::MPI.Comm) :: Vector{Hamiltonian}
	if root == MPI.Comm_rank(comm)
			@assert Hs != nothing
			len = length(Hs)
			n_procs = MPI.Comm_size(comm)
			reqs = Vector{MPI.Request}(undef, n_procs - 1)
			for i in 1:n_procs-1
					reqs[i] = MPI.isend(Hs[p_range(len, n_procs, i+1)], i, 0, comm)
			end
			qo_self = Hs[p_range(len, n_procs, 1)]
			MPI.Waitall!(reqs)
			return qo_self
	else
			qo, status = MPI.recv(root, 0, comm)
			return qo
	end
end