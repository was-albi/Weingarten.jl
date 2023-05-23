module Weingarten

# actual export
export WgU, WgO, compute_unitary_expectation, compute_orthogonal_expectation
# for debugging purposes
export Id_Hn, Gelfand_product, Gelfand_convolution, zonal_spherical_function

function C_lambda(l, z)
    out = 1
    for i in 1:length(l)
        for j in 1:l[i]
            out = out*(z + j -i)
        end
    end
    return out
end

function WgU(sigma, d)
    out = 0
    sigmatype = permtype(sigma)
    n = sum(sigmatype)
    lambdas = Generic.partitions(n)
    for l in 1:length(lambdas)
        c = C_lambda(lambdas[l],d)
        if abs(c) > 10^(-12)
            chi = character(lambdas[l])
            out = out + chi(Perm(n))* chi(sigma)/c
        end
    end
    return out/factorial(n)
end

function multi_delta(sigma, multi_ii, multi_jj)
    out_ = 1
    for it in 1:length(multi_ii)
        if multi_ii[getindex(sigma, it)] != multi_jj[it]
            out_ = 0
        end
    end

    return out_
end

function compute_unitary_expectation(d, multi_i, multi_j, multi_k, multi_l)
    # compute expectation value of product of elements in U (unitary) according to multiindices i,j,k,l
    out = 0
    n = length(multi_i)
    if n != length(multi_j) || n != length(multi_j) || n != length(multi_j)
        println("Error: multi indices have different length")
        return "Error"
    end

    sigmas = SymmetricGroup(n)
    taus = SymmetricGroup(n)

    for s in sigmas
        for t in taus

            out = out + multi_delta(s,multi_i,multi_k)*multi_delta(t,multi_j,multi_l)*WgU(inv(s)*t,d)
            #out = out + multi_delta(s,multi_i,multi_k)*multi_delta(t,multi_j,multi_l)*WgU(inv(t)*s,d)

        end
    end

    return out

end

function check_M2n(sigma, n)
    flag = 1
    for i in 1:(n-1)
        if getindex(sigma, 2*i-1) >= getindex(sigma, 2*i)
            flag = 0
        end
        if getindex(sigma, 2*i-1) >= getindex(sigma, 2*i+1)
            flag = 0
        end
    end
    if getindex(sigma, 2*n-1) >= getindex(sigma, 2*n)
        flag = 0
    end

    return flag
end


function Id_Hn(sigma, n)
    # Q: can I set t as a const instead of re-computing it every time? It depends on n tho
    t = Perm([2;1; collect(3:2*n)])
    for it in 3:2:(2*n-1)
        s = Perm([collect(1:it-1);it+1;it;collect(it+2:2*n)])
        t = t*s
    end
    flag = 0
    if sigma*t == t*sigma
        flag = 1
    end
    return flag
end


function Gelfand_product(f1,f2,n)
    # return f1 star f2
    S2n = SymmetricGroup(2*n)
    #M2n = []
    #for sigma in S2n
    #    if check_M2n(sigma, n) == 1
    #        push!(M2n, sigma)
    #    end
    #end
    function f1_star_f2(sigma)
        S2n = SymmetricGroup(2*n)
        out = 0
        for tau in S2n
            if check_M2n(tau,n) == 1
                out += f1(tau)*f2(inv(tau)*sigma)
            end
        end
        return out
    end

    return f1_star_f2
end

function Gelfand_convolution(f1,f2,n)

    function f1_asterix_f2(sigma)
        f1_star_f2 = Gelfand_product(f1,f2,n)
        return f1_star_f2(sigma)*(2^n*factorial(n))
    end

    return f1_asterix_f2
end


function zonal_spherical_function(sigma, lambda_)
    n = sum(lambda_)
    chi = character(Partition(2*lambda_))
    S2n = SymmetricGroup(2*n)
    out = 0
    for s in S2n
        out += Id_Hn(s, n)*chi(sigma*s)
    end
    return out/(2^n*factorial(n))
end


function D_lambda(l, z)
    out = 1
    for i in 1:length(l)
        for j in 1:l[i]
            out = out*(z + 2*j - i - 1)
        end
    end
    return out
end


function WgO(sigma, d)
    out = 0
    sigmatype = permtype(sigma)
    n = div(sum(sigmatype),2)
    lambdas = Generic.partitions(n)
    for l in 1:length(lambdas)
        D = D_lambda(lambdas[l],d)
        if abs(D) > 10^(-12)
            chi = character(Partition(2*lambdas[l]))
            out = out + chi(Perm(2*n))* zonal_spherical_function(sigma, lambdas[l])/D
        end
    end
    return out* (2^n*factorial(n)/factorial(2*n))
end


function multi_Big_delta(sigma, multi_ii)
    out_ = 1
    n = div(length(multi_ii),2)
    for it in 1:n
        if multi_ii[getindex(sigma, 2*it-1)] != multi_ii[getindex(sigma, 2*it)]
            out_ = 0
        end
    end

    return out_
end


function compute_orthogonal_expectation(d, multi_i, multi_j)
    # compute expectation value of product of elements in U (unitary) according to multiindices i,j
    out = 0
    nn = length(multi_i)
    if nn != length(multi_j)
        println("Error: multi indices have different lengths")
        return "Error"
    end
    if mod(nn,2) != 0
        return 0
    end

    n = div(nn, 2)

    S2n = SymmetricGroup(nn)
    M2n = []
    for sigma in S2n
        if check_M2n(sigma, n) == 1
            push!(M2n, sigma)
        end
    end

    for s in M2n
        for t in M2n

            out = out + multi_Big_delta(s,multi_i)*multi_Big_delta(t,multi_j)*WgO(inv(s)*t,d)

        end
    end

    return out

end

end
