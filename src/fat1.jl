export lu_mark

"""
    lu_mark(A)
    Fatoração LU de matriz esparsa usando algoritmo de Markowitz
"""
###
function lu_mark(A,debug=false)
    # L: matriz triangular inferior unitária
    # U: matriz triangular superior
    # p: vetor de permutação de linhas
    # q: vetor de permutação de colunas

    m,n = size(A)

    L = zeros(m,n)
    U = zeros(m,n)
    B = Array{Float64}(A[:,:])
    p = collect(1:m)
    q = collect(1:n)

    for k = 1:n-1
        debug && println("k:",k)
        # Heurística de busca de pivô
        r = Vector{Int64}() # Quantidade de elementos não zerados na linha
        c = Vector{Int64}() # Quantidade de elementos não zerados na coluna
        # Varrendo linha a linha da submatriz
        for i = k:m
            rr = 0
            # Varrendo coluna a coluna da submatriz
            for j = k:n
                # Contando elementos não zerados
                if( B[i,j] != 0 )
                    rr += 1
                end
            end
            # Colocando na lista
            append!(r,rr)
        end
        # Varrendo coluna a coluna da submatriz
        for j = k:n
            cc = 0
            # Varrendo linha a linha da submatriz
            for i = k:m
                # Contando elementos não zerados
                if( B[i,j] != 0 )
                    cc += 1
                end
            end
            # Colocando na lista
            append!(c,cc)
        end
        debug && println("B:")
        debug && display(B)
        debug && println("r:",r)
        debug && println("c:",c)
        # Agora vamos encontrar o mínimo (r[i]-1)*(c[j]-1)
        r_escolhido = 1
        c_escolhido = 1
        menor = (r[1]-1)*(c[1]-1)
        aha = false
        for j = 1:length(c)
            for i = 1:length(r)
                #println("[",k+i-1,",",k+j-1,"]")
                if( ( B[k-1+i,k-1+j] != 0 ) && ( (r[i]-1)*(c[j]-1) < menor ))
                    r_escolhido = i
                    c_escolhido = j
                    menor = (r[i]*c[j])
                    if( r[i] == 1 || c[j] == 1 )
                        aha = true
                        #println("aha")
                        break
                    end
                end
            end
            if( aha )
                break
            end
        end
        debug && println("r_e:",r_escolhido,";c_e:",c_escolhido)
        if( r_escolhido != 1 )
            # Fazendo a troca de linhas
            for j = 1:n
                B[k,j],B[k+r_escolhido-1,j] = B[k+r_escolhido-1,j],B[k,j]
            end
            
            for j = 1:k-1
                L[k,j],L[k+r_escolhido-1,j] = L[k+r_escolhido-1,j],L[k,j]
            end

            # Armazenando a informação de troca de linhas
            p[k],p[k+r_escolhido-1] = p[k+r_escolhido-1],p[k]
        end
        if( c_escolhido != 1 )
            # Fazendo a troca de colunas
            for i = 1:n
                B[i,k],B[i,k+c_escolhido-1] = B[i,k+c_escolhido-1],B[i,k]
            end
            # Armazenando a informação de troca de colunas
            q[k],q[k+c_escolhido-1] = q[k+c_escolhido-1],q[k]
        end
        debug && println("p':", p)
        debug && println("q':", q)
        debug && println("B':")
        debug && display(B)
        # Eliminação Gaussiana considerando que os pivôs estão na diagonal principal
        L[k,k] = 1
        for i = k+1:m
            if( B[i,k] != 0)
                mik = B[i,k] / B[k,k]
                L[i,k] = mik
                for j = k:n
                    B[i,j] = B[i,j] - B[k,j]*mik 
                end
            end
        end
        # Populando elementos da matriz U
        for i = 1:k
            U[i,k] = B[i,k]
        end
        debug && println("L':")
        debug && display(L)
        debug && println("U':")
        debug && display(U)
        
    end
    debug && println("B'':")
    debug && display(B)
    L[n,n] = 1
    for i = 1:m
        U[i,n] = B[i,n]
    end
    return L,U,p,q
end


