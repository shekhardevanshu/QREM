using Printf
using DelimitedFiles, Statistics
using FLoops
using LinearAlgebra
# using KernelDensity
# BLAS.set_num_threads(1)

function binaryRead(fname)
    data = Array{Float64}(undef, filesize(fname)÷8)
    read!(fname, data)
    # data = ntoh.(data)
    # data_f = filter(x -> abs(x) > tol, data)
    return data
end

function getSpacingRatio(eVal) #Phys. Rev. B 75, 155111

    nconv = length(eVal)
    spr_s = []
    
    for n in 2:nconv-1
        dn = eVal[n+1] - eVal[n]
        dn_1 = eVal[n] - eVal[n-1]
        
        if dn != 0 && dn_1 != 0
            push!(spr_s, ( min(dn, dn_1) / max(dn, dn_1) ))
        # else
        #     print("eVal[n+1] $(eVal[n+1]), eVal[n] $(eVal[n])")
        #     flush(stdout)
        end

    end

    return mean(spr_s)

end

function main(args)
    
    # barr = sort(vcat(collect(0.0:1:5), collect(10:5:20), collect(30:20:100),  collect(150:50:256), [120, 180, 500], collect(7:2:14), collect(0.1:0.2:1))) # 10, 12
    # barr = sort(vcat([0.1], collect(1:1:5), collect(10:5:20), collect(30:20:100), collect(35:20:100), collect(150:50:256), [120, 180], collect(500:200:1500), collect(7:2:14))) # 14
    barr = sort(vcat([0.1], collect(1:1:5), collect(10:5:20), collect(30:20:100), collect(25:10:50), collect(150:50:256), [120, 180, 350, 450], collect(500:200:1500), collect(7:2:14))) # 16
    # barr = sort(vcat(collect(0.1:0.2:10), collect(10:1:30), collect(35:5:100))) # 12
    bname = map(x -> @sprintf("%0.1f", x), barr)
   
    l = parse(Int, args[1])
    s = parse(Int, args[3])
    n_eval = parse(Int, args[2])
    # a = parse(Float64, args[4])
    # c = parse(Float64, args[5])
    # c = @sprintf("%0.1f", c)
    # a = @sprintf("%0.2f", a)
    e = "0.0"

    spr_h = Matrix{Float64}(undef, length(barr), s)
    lsp_h = Matrix{Float64}(undef, length(barr), s*(n_eval-1))
    
    spr_s = zeros(s)
    lsp_s = zeros(s*(n_eval-1))
    
    for (j, g) in enumerate(bname)

        for m in 1:s
            
		eVal = binaryRead("./dataBinQREM_L=$(l)/eigval_L=$(l)_b=$(g)_itr_$(m)_e=$e.dat")
		# eVal = binaryRead("./dataBinQREM_L=$(l)/eigval_L=$(l)_b=$(g)_a=$(a)_g=$(c)_itr_$(m)_e=$e.dat")
		eVal2 = sort!(eVal[1:n_eval])
            
                lsp_s[(n_eval-1)*(m-1)+1:(n_eval-1)*m] = diff(eVal2)
	        spr_s[m] = getSpacingRatio(eVal2)
            
        end

	print("completed g = $g, mlsp = $(mean(lsp_s))\n")
	print("completed g = $g, mspr = $(mean(spr_s))\n")
        flush(stdout)

        spr_h[j, :] = spr_s
        lsp_h[j,:] = lsp_s
        
    end

    writedlm("lsp_qrem_pe,nev=$n_eval,L=$l,e=$e.txt", [lsp_h barr])
    writedlm("spr_qrem_pe,nev=$n_eval,L=$l,e=$e.txt", [spr_h barr])
end

main(ARGS)
