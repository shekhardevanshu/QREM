using Printf
using DelimitedFiles, Statistics
using QuantumInformation
using FLoops
using LinearAlgebra

BLAS.set_num_threads(1)

function binaryRead(fname)
    data = Array{Float64}(undef, filesize(fname)÷8)
    read!(fname, data)
    return data
end

function get_schmidt(U, d)
    C = zeros(ComplexF64, d, d)
    rho = similar(C)
    C = reshape(U, d, d)
    mul!(rho, C, C')
    return vonneumann_entropy(rho)/log(2), tr(rho^2)
end

function get_entropy(U, d, nconv, N)
    vn_phi = Vector{Float64}(undef, nconv)
    s2_phi = Vector{Float64}(undef, nconv)
    psi = zeros(ComplexF64, N)

    for n in 1:nconv

        psi .= U[(n-1)*N+1:n*N]
	vn_phi[n], s2_phi[n] = get_schmidt(psi, d)

    end

    vn_phi, s2_phi
end

function main(arbs)
    L = parse(Int, arbs[1])
    nconv = parse(Int, arbs[2])
    s = parse(Int, arbs[3])
    # a = parse(Float64, arbs[4])
    # g = parse(Float64, arbs[5])
    # a = @sprintf("%0.2f", a)
    # g = @sprintf("%0.1f", g)

    e = "0.0"

    # barr = sort(vcat(collect(0.0:1:5), collect(10:5:20), collect(30:20:100),  collect(150:50:256), [120, 180, 500], collect(7:2:14), collect(0.1:0.2:1))) # 10, 12
    # barr = sort(vcat([0.1], collect(1:1:5), collect(10:5:20), collect(30:20:100), collect(35:20:100), collect(150:50:256), [120, 180], collect(500:200:1500), collect(7:2:14))) # 14
    barr = sort(vcat([0.1], collect(1:1:5), collect(10:5:20), collect(30:20:100), collect(25:10:50), collect(150:50:256), [120, 180, 350, 450], collect(500:200:1500), collect(7:2:14))) # 16
    # barr = sort(vcat(collect(0.1:0.2:10), collect(10:1:30), collect(35:5:100))) # 12
    bname = map(x -> @sprintf("%0.1f", x), barr)
    N = 2^L
    d = 2^(div(L, 2))
    
    vn_s = Array{Float64}(undef, s*nconv)
    vn_h = Matrix{Float64}(undef, length(barr), length(vn_s))
    s2_s = Array{Float64}(undef, s*nconv)
    s2_h = Matrix{Float64}(undef, length(barr), length(vn_s))

    for (j, b) in enumerate(bname)

        @floop for m in 1:s

	    eVec = binaryRead("./dataBinQREM_L=$L/eigvec_L=$(L)_b=$(b)_itr_$(m)_e=$e.dat")
	    # eVec = binaryRead("./dataBinQREM_L=$L/eigvec_L=$(L)_b=$(b)_a=$(a)_g=$(g)_itr_$(m)_e=$e.dat")
            vn_s[(m-1)*nconv+1:nconv*m], s2_s[(m-1)*nconv+1:nconv*m] = get_entropy(eVec, d, nconv, N)

          
        end

	print("completed b = $b, mean = $(mean(s2_s)), var = $(var(s2_s))\n")
        flush(stdout)
        vn_h[j,:] = vn_s
	s2_h[j,:] = s2_s
                
    end

    # writedlm("vn_qrem_pe,L=$L,nev=$nconv,a=$a,g=$g,e=$e.txt", vn_h)
    # writedlm("s2_qrem_pe,L=$L,nev=$nconv,a=$a,g=$g,e=$e.txt", s2_h)
   
    writedlm("vn_qrem_pe,L=$L,nev=$nconv,e=$e.txt", [vn_h barr])
    writedlm("s2_qrem_pe,L=$L,nev=$nconv,e=$e.txt", [s2_h barr])
end

main(ARGS)
