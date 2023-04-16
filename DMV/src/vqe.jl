

function get_expectation(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observable::Hamiltonian, trunc::MPSTruncation)
    temp = normalize_mps(apply(reset_parameters_wrapper(apply_map(θs, map=map), x=circuit), ψ, trunc=trunc))
    #println(QuantumSpins.bond_dimension(temp))
    real(expectation(observable, temp, trunc=trunc))
end

function get_expectation(ψ::MPS, observable::Hamiltonian, trunc::MPSTruncation)
    real(expectation(observable, ψ, trunc=trunc))
end


function get_expectation_and_gradient_no_timing(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observable::Hamiltonian, trunc::MPSTruncation)
    if length(observable.data) == 0
        return zeros(Float64, length(θs) + 1), zeros(Float64, 2)
    end
    obj(x) = get_expectation(x, map, circuit, ψ, observable, trunc)
    value, back = pullback(obj, θs)
    pack = back(1.)[1]
    pushfirst!(pack, value)
    return pack, zeros(Float64, 2)
end

function get_expectation_and_gradient_timing(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observable::Hamiltonian, trunc::MPSTruncation)
    if length(observable.data) == 0
        return zeros(Float64, length(θs) + 1), zeros(Float64, 2)
    end
    
    obj(x) = get_expectation(x, map, circuit, ψ, observable, trunc)

    t1 = time_ns()
    value, back = pullback(obj, θs)
    t2 = time_ns()
    pack = back(1.)[1]
    t3 = time_ns()

    pushfirst!(pack, value)
    return pack, [(t2-t1)/1e9, (t3-t2)/1e9]
end

get_expectation_and_gradient(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observable::Hamiltonian, trunc::MPSTruncation
                                ; timing=false) = timing ? 
                                    get_expectation_and_gradient_timing(θs, map, circuit, ψ, observable, trunc) : 
                                    get_expectation_and_gradient_no_timing(θs, map, circuit, ψ, observable, trunc)



function get_expectation(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observables::Vector{Hamiltonian}, trunc::MPSTruncation)
    v = 0.0
    for i in 1:length(observables)
        v += get_expectation(θs, map, circuit, ψ, observables[i], trunc)
    end
    v
end

function get_expectation(ψ::MPS, observables::Vector{Hamiltonian}, trunc::MPSTruncation)
    v = 0.0
    for i in 1:length(observables)
        v += get_expectation(ψ, observables[i], trunc)
    end
    v
end

function get_expectation_and_gradient_no_timing(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observables::Vector{Hamiltonian}, trunc::MPSTruncation)
    pack = zeros(length(θs) + 1)
    for i in 1:length(observables)
        if length(observables[i].data) > 0
            obj(x) = get_expectation(x, map, circuit, ψ, observables[i], trunc)
            value, back = pullback(obj, θs)
            pack[2:end] += back(1.)[1]
            pack[1] += value
        end
    end
    return pack, [0.0, 0.0]
end

function get_expectation_and_gradient_timing(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observables::Vector{Hamiltonian}, trunc::MPSTruncation)
    pack = zeros(length(θs) + 1)
    total_t1 = 0
    total_t2 = 0

    for i in 1:length(observables)
        if length(observables[i].data) > 0
            obj(x) = get_expectation(x, map, circuit, ψ, observables[i], trunc)

            t1 = time_ns()
            value, back = pullback(obj, θs)
            t2 = time_ns()
            pack[2:end] += back(1.)[1]
            t3 = time_ns()

            pack[1] += value
            total_t1 += (t2 - t1)
            total_t2 += (t3 - t2)
        end
    end
    return pack, [total_t1/1e9, total_t2/1e9]
end

get_expectation_and_gradient(θs::Vector, map::Vector, circuit::QCircuit, ψ::MPS, observables::Vector{Hamiltonian}, trunc::MPSTruncation
                                ; timing=false) = timing ? 
                                    get_expectation_and_gradient_timing(θs, map, circuit, ψ, observables, trunc) : 
                                    get_expectation_and_gradient_no_timing(θs, map, circuit, ψ, observables, trunc)