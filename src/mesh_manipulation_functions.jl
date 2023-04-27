using SparseArrays

function From_3D_to_1D(i, j, k, M, N)
    pos = ((k-1) * M * N) + ((j-1) * M) + i;
    return pos
end

function bin_search(num, A)

    index = 1
    n = size(A)[1]
    left = 1
    right = n

    while left <= right
        mid = ceil(Int, (left + right) / 2)
        #mid = Int64((left + right) / 2)

        if A[mid] == num
            index = mid
            break
        else
            if A[mid] > num
                right = mid - 1
            else
                left = mid + 1
            end
        end
    end
    return index
end

function create_volumes_mapping_and_centers(matrice,Nx,Ny,Nz,num_centri,sx,sy,sz,min_v)


    println("----",(Nx * Ny * Nz))
 
    mapping = zeros(Int64, Nx * Ny * Nz)
    centri_vox = zeros(Float64, num_centri, 3)
    id_mat = zeros(Int64, num_centri)

    num_grids = size(matrice)[1]
    num_ele=1
    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3]==true)
                        mapping[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]=num_ele
                        centri_vox[num_ele, 1] = min_v[1] + (sx * (cont - 1))  + (sx / 2.0)
                        centri_vox[num_ele, 2] = min_v[2] + (sy * (cont2 - 1)) + (sy / 2.0)
                        centri_vox[num_ele, 3] = min_v[3] + (sz * (cont3 - 1)) + (sz / 2.0)
                        id_mat[num_ele] = k
                        num_ele = num_ele + 1
                        break
                    end
                end
            end
        end
    end


    return num_ele,mapping,centri_vox,id_mat
end

function create_nodes_ref(matrice, Nx,Ny,Nz, num_centri, external_g, m_volumes)

    num_grids = size(matrice)[1]
    nodes = zeros(Int64, num_centri)

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        c1 = 1 + 2 * (cont - 1) + 1
                        c2 = 1 + 2 * (cont2 - 1) + 1
                        c3 = 1 + 2 * (cont3 - 1) + 1
                        f1 = external_g[k, 1, cont, cont2, cont3]
                        f2 = external_g[k, 2, cont, cont2, cont3]
                        f3 = external_g[k, 3, cont, cont2, cont3]
                        f4 = external_g[k, 4, cont, cont2, cont3]
                        f5 = external_g[k, 5, cont, cont2, cont3]
                        f6 = external_g[k, 6, cont, cont2, cont3]
                        is_f1 = f1
                        is_f2 = f2
                        if (f1==1 && f2==1)
                            is_f1 = 0
                            is_f2 = 0
                        end

                        is_f3 = f3
                        is_f4 = f4
                        if (f3==1 && f4==1)
                            is_f3 = 0
                            is_f4 = 0
                        end

                        is_f5 = f5
                        is_f6 = f6
                        if (f5 == 1 && f6 == 1)
                            is_f5 = 0
                            is_f6 = 0
                        end

                        if (is_f1==1 || is_f2==1 || is_f3==1 || is_f4==1 || is_f5==1 || is_f6==1) 
                            if (is_f1==1 && is_f2==0 && is_f3==0 && is_f4==0 && is_f5==0 && is_f6==0) 
                                nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                    From_3D_to_1D(c1, c2 - 1, c3, 3 * Nx, 3 * Ny)
                            else
                                if (is_f2 == 1 && is_f1 == 0 && is_f3 == 0 && is_f4 == 0 && is_f5 == 0 && is_f6 == 0)
                                    nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                        From_3D_to_1D(c1, c2 + 1, c3, 3 * Nx, 3 * Ny)
                                else
                                    if (is_f3 == 1 && is_f1 == 0 && is_f2 == 0 && is_f4 == 0 && is_f5 == 0 && is_f6 == 0)
                                        nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                            From_3D_to_1D(c1-1, c2 , c3, 3 * Nx, 3 * Ny)
                                    else
                                        if (is_f4 == 1 && is_f1 == 0 && is_f2 == 0 && is_f3 == 0 && is_f5 == 0 && is_f6 == 0)
                                            nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                From_3D_to_1D(c1 + 1, c2, c3, 3 * Nx, 3 * Ny)
                                        else
                                            if (is_f5 == 1 && is_f1 == 0 && is_f2 == 0 && is_f3 == 0 && is_f4 == 0 && is_f6 == 0)
                                                nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                    From_3D_to_1D(c1, c2, c3-1, 3 * Nx, 3 * Ny)
                                            else
                                                if (is_f6 == 1 && is_f1 == 0 && is_f2 == 0 && is_f3 == 0 && is_f4 == 0 && is_f5 == 0)
                                                    nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                        From_3D_to_1D(c1, c2, c3 + 1, 3 * Nx, 3 * Ny)
                                                else
                                                    if (is_f1 == 1 && is_f3 == 1 && is_f2 == 0 && is_f4 == 0 && is_f5 == 0 && is_f6 == 0)
                                                        nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                            From_3D_to_1D(c1-1, c2-1, c3, 3 * Nx, 3 * Ny)
                                                    else
                                                        if (is_f1 == 1 && is_f4 == 1 && is_f2 == 0 && is_f3 == 0 && is_f5 == 0 && is_f6 == 0)
                                                            nodes[m_volumes[
                                                                From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                                From_3D_to_1D(c1 + 1, c2 - 1, c3, 3 * Nx, 3 * Ny)
                                                        else
                                                            if (is_f1 == 1 && is_f5 == 1 && is_f2 == 0 && is_f3 == 0 && is_f4 == 0 && is_f6 == 0)
                                                                nodes[m_volumes[
                                                                    From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                                    From_3D_to_1D(c1, c2 - 1, c3 - 1, 3 * Nx, 3 * Ny)
                                                            else
                                                                if (is_f1 == 1 && is_f6 == 1 && is_f2 == 0 && is_f3 == 0 && is_f4 == 0 && is_f5 == 0)
                                                                    nodes[m_volumes[
                                                                        From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]] = 
                                                                        From_3D_to_1D(c1, c2 - 1, c3 + 1, 3 * Nx,3 * Ny)
                                                                else
                                                                    if (is_f1 == 1 && is_f3 == 1 && is_f5 == 1 && is_f2 == 0 && is_f4 == 0 && is_f6 == 0)
                                                                        nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx,Ny)]] = 
                                                                            From_3D_to_1D(c1-1, c2 - 1, c3 - 1, 3 * Nx, 3 * Ny)
                                                                    else
                                                                        if (is_f1 == 1 && is_f3 == 1 && is_f6 == 1 && is_f2 == 0 && is_f4 == 0 && is_f5 == 0)
                                                                            nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3, Nx,Ny)]] = 
                                                                                From_3D_to_1D(c1 - 1, c2 - 1, c3 + 1,3 * Nx, 3 * Ny)
                                                                        else
                                                                            if (is_f1 == 1 && is_f4 == 1 && is_f5 == 1 && is_f2 == 0 && is_f3 == 0 && is_f6 == 0)
                                                                                nodes[m_volumes[From_3D_to_1D(cont, cont2, cont3,Nx, Ny)]] = 
                                                                                    From_3D_to_1D(c1 + 1, c2 - 1,c3 - 1, 3 * Nx,3 * Ny)
                                                                            else
                                                                                if (is_f1 == 1 && is_f4 == 1 && is_f6 == 1 && is_f2 == 0 && is_f3 == 0 && is_f5 == 0)
                                                                                    nodes[m_volumes[From_3D_to_1D(cont, cont2,cont3, Nx, Ny)]] = 
                                                                                        From_3D_to_1D(c1 + 1, c2 - 1,c3 + 1, 3 * Nx,3 * Ny)
                                                                                else
                                                                                    if (is_f2 == 1 && is_f3 == 1 && is_f1 == 0 && is_f4 == 0 && is_f5 == 0 && is_f6 == 0)
                                                                                        nodes[m_volumes[From_3D_to_1D(cont, cont2,cont3, Nx,Ny)]] = 
                                                                                            From_3D_to_1D(c1 - 1,c2 + 1,c3,3 * Nx,3 * Ny)
                                                                                    else
                                                                                        if (is_f2 == 1 && is_f4 == 1 && is_f1 == 0 && is_f3 == 0 && is_f5 == 0 && is_f6 == 0)
                                                                                            nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3, Nx,Ny)]] = 
                                                                                                From_3D_to_1D(c1 + 1,c2 + 1,c3,3 * Nx,3 * Ny)
                                                                                        else
                                                                                            if (is_f2 == 1 && is_f5 == 1 && is_f1 == 0 && is_f3 == 0 && is_f4 == 0 && is_f6 == 0)
                                                                                                nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                    From_3D_to_1D(c1, c2 + 1,c3 - 1, 3 * Nx,3 * Ny)
                                                                                            else
                                                                                                if (is_f2 == 1 && is_f6 == 1 && is_f1 == 0 && is_f3 == 0 && is_f4 == 0 && is_f5 == 0)
                                                                                                    nodes[m_volumes[From_3D_to_1D(cont, cont2,cont3, Nx,Ny)]] = 
                                                                                                        From_3D_to_1D(c1, c2 + 1,c3 + 1,3 * Nx,3 * Ny)
                                                                                                else
                                                                                                    if (is_f2 == 1 && is_f3 == 1 && is_f5 == 1 && is_f1 == 0 && is_f4 == 0 && is_f6 == 0)
                                                                                                        nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                            From_3D_to_1D(c1 - 1,c2 + 1,c3 - 1,3 * Nx,3 * Ny)
                                                                                                    else
                                                                                                        if (is_f2 == 1 && is_f3 == 1 && is_f6 == 1 && is_f1 == 0 && is_f4 == 0 && is_f5 == 0)
                                                                                                            nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                                From_3D_to_1D(c1 - 1,c2 + 1,c3 + 1, 3 * Nx, 3 * Ny)
                                                                                                        else
                                                                                                            if (is_f2 == 1 && is_f4 == 1 && is_f5 == 1 && is_f1 == 0 && is_f3 == 0 && is_f6 == 0)
                                                                                                                nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                                    From_3D_to_1D(c1 + 1,c2 + 1,c3 - 1,3 * Nx,3 * Ny)
                                                                                                            else
                                                                                                                if (is_f2 == 1 && is_f4 == 1 && is_f6 == 1 && is_f1 == 0 && is_f3 == 0 && is_f5 == 0)
                                                                                                                    nodes[m_volumes[From_3D_to_1D( cont,cont2,cont3,Nx, Ny)]] = 
                                                                                                                        From_3D_to_1D(c1 + 1,c2 + 1,c3 + 1,3 * Nx, 3 * Ny)
                                                                                                                else
                                                                                                                    if (is_f3 == 1 && is_f5 == 1 && is_f1 == 0 && is_f2 == 0 && is_f4 == 0 && is_f6 == 0)
                                                                                                                        nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                                            From_3D_to_1D(c1 - 1,c2 ,c3 - 1,3 * Nx,3 * Ny)
                                                                                                                    else
                                                                                                                        if (is_f3 == 1 && is_f6 == 1 && is_f1 == 0 && is_f2 == 0 && is_f4 == 0 && is_f5 == 0)
                                                                                                                            nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                                                From_3D_to_1D(c1 - 1,c2, c3 + 1,3 * Nx,3 * Ny)
                                                                                                                        else
                                                                                                                            if (is_f4 == 1 && is_f5 == 1 && is_f1 == 0 && is_f2 == 0 && is_f3 == 0 && is_f6 == 0)
                                                                                                                                nodes[m_volumes[From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                                                    From_3D_to_1D(c1 + 1,c2,c3 - 1,3 * Nx,3 * Ny)
                                                                                                                            else
                                                                                                                                if (is_f4 == 1 && is_f6 == 1 && is_f1 == 0 && is_f2 == 0 && is_f3 == 0 && is_f5 == 0)
                                                                                                                                    nodes[m_volumes[ From_3D_to_1D(cont,cont2,cont3,Nx,Ny)]] = 
                                                                                                                                    From_3D_to_1D(c1 + 1,c2,c3 + 1,3 * Nx,3 * Ny)
                                                                                                                                end
                                                                                                                            end
                                                                                                                        end
                                                                                                                    end
                                                                                                                end
                                                                                                            end
                                                                                                        end
                                                                                                    end
                                                                                                end
                                                                                            end
                                                                                        end
                                                                                    end
                                                                                end
                                                                            end
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            nodes[m_volumes[From_3D_to_1D( cont,cont2,cont3,Nx, Ny)]] = From_3D_to_1D(c1,c2,c3,3 * Nx, 3 * Ny)
                        end
                        break
                    end
                end
            end
        end
    end
            
    nodes_red = sort(unique(nodes))

    return nodes_red,nodes
end

function slicematrix(A::AbstractMatrix{T}) where T
    m, n = size(A)
    B = Vector{T}[Vector{T}(undef, n) for _ in 1:m]
    for i in 1:m
        B[i] .= A[i, :]
    end
    return B
end

function distfcm(center, data)
    function sum_el(array)
        sum = 0
        for el in array
            sum = sum + el
        end
        return sum
    end

    out = slicematrix(sum(sum_el, (data .- (ones(size(data)[1],1) .* center)[1][1]).^2, dims=2).^0.5)
    return out
end

function nodes_find_rev(Nodes_inp_coord, nodi_centri, node_to_skip)

    indici = sortperm(distfcm(Nodes_inp_coord, nodi_centri))

    if (indici[1]!=node_to_skip)
        node = indici[1]
    else
        node = indici[2]
    end

    return node
end

function find_nodes_port(nodi_centri, port_start, port_end, nodi, nodi_red)

    N = size(port_start)[1]
    port_voxels = zeros(Int64, N, 2)
    port_nodes = zeros(Int64, N, 2)

    for cont in range(1, stop=N)
        port_voxels[cont, 1] = nodes_find_rev(port_start[cont,:], nodi_centri, -1)
        port_voxels[cont, 2] = nodes_find_rev(port_end[cont,:],nodi_centri, port_voxels[cont,1])
        port_nodes[cont, 1] = bin_search(nodi[port_voxels[cont, 1]],nodi_red)
        port_nodes[cont, 2] = bin_search(nodi[port_voxels[cont, 2]],nodi_red)
    end



    return port_voxels, port_nodes
end

function create_external_grids(matrice,Nx,Ny,Nz)
    num_grids = size(matrice)[1]
    
    OUTPUTgrids = zeros(Int8 ,num_grids ,6,Nx ,Ny ,Nz)


    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][1][cont2][cont3]==1)
                    OUTPUTgrids[k,3,1,cont2,cont3]=1
                end
            end
        end
    end

    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][Nx][cont2][cont3]==1)
                    OUTPUTgrids[k,4,Nx,cont2,cont3]=1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][cont][1][cont3]==1)
                    OUTPUTgrids[k,1,cont,1,cont3] = 1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if (matrice[k][cont][Ny][cont3] == 1)
                    OUTPUTgrids[k,2,cont,Ny,cont3] = 1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if (matrice[k][cont][cont2][1]==1)
                    OUTPUTgrids[k,5,cont,cont2,1] = 1
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if(matrice[k][cont][cont2][Nz]==1)
                    OUTPUTgrids[k,6,cont,cont2,Nz] = 1
                end
            end
        end
    end

    for cont in range(2,stop=Nx-1)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        if (matrice[k][cont - 1][cont2][cont3] == 0)
                            OUTPUTgrids[k,3,cont,cont2,cont3] = 1
                        end
                        if (matrice[k][cont + 1][cont2][cont3] == 0)
                            OUTPUTgrids[k,4,cont,cont2,cont3] = 1
                        end
                    end
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(2, stop=Ny-1)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        if (matrice[k][cont][cont2 - 1][cont3] == 0)
                            OUTPUTgrids[k,1,cont,cont2,cont3] = 1
                        end
                        if (matrice[k][cont][cont2 + 1][cont3] == 0)
                            OUTPUTgrids[k,2,cont,cont2,cont3] = 1
                        end
                    end
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(2, stop=Nz-1)
                for k in range(1, stop=num_grids)
                    if (matrice[k][cont][cont2][cont3] == 1)
                        if (matrice[k][cont][cont2][cont3 - 1] == 0)
                            OUTPUTgrids[k,5,cont,cont2,cont3] = 1
                        end
                        if (matrice[k][cont][cont2][cont3 + 1] == 0)
                            OUTPUTgrids[k,6,cont,cont2,cont3] = 1
                        end
                    end
                end
            end
        end
    end

    return OUTPUTgrids
end

function create_mapping_Ax(matrice,Nx,Ny,Nz)
    num_grids = size(matrice)[1]

    N_max = ((Nx - 1) * Ny * Nz)
    mapping = zeros(Int64, N_max)

    num_ele = 0

    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for cont in range(1, stop=Nx-1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont+1][cont2][cont3]==1))
                        kkey = From_3D_to_1D(cont, cont2, cont3, Nx - 1, Ny)
                        if mapping[kkey] == 0
                            mapping[kkey] = num_ele
                            num_ele = num_ele + 1
                        end
                        break
                    end
                end
            end
        end
    end

    return mapping,num_ele
end

function create_mapping_Ay(matrice,Nx,Ny,Nz)
    num_grids = size(matrice)[1]

    N_max = ( Nx * (Ny - 1) * Nz)
    mapping = zeros(Int64 ,N_max)

    num_ele = 0

    for cont3 in range(1, stop=Nz)
        for cont in range(1, stop=Nx)
            for cont2 in range(1, stop=Ny - 1)
                for k in range(1, stop=num_grids)
                    if((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont][cont2+1][cont3]==1))
                        kkey = From_3D_to_1D(cont, cont2, cont3, Nx, Ny - 1)
                        if mapping[kkey] == 0
                           mapping[kkey] = num_ele
                           num_ele = num_ele + 1
                        end
                        break
                    end
                end
            end
        end
    end

    return mapping,num_ele
end

function create_mapping_Az(matrice,Nx,Ny,Nz)
    num_grids = size(matrice)[1]
    N_max = ( Nx * Ny * (Nz - 1) )
    mapping = zeros(Int64, N_max)

    num_ele = 0

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz - 1)
                for k in range(1, stop=num_grids)
                    if((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont][cont2][cont3+1]==1))
                        kkey = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        if mapping[kkey] == 0
                           mapping[kkey] = num_ele
                           num_ele = num_ele + 1
                        end
                        break
                    end
                end
            end
        end
    end

    return mapping,num_ele
end

function create_A_mats_volInd(matrice,Nx,Ny,Nz,mapping_Vox,mapAx, NAx, mapAy, NAy, mapAz, NAz, sx,sy,sz,min_v,nodi,nodi_red)
 
    num_grids = size(matrice)[1]
    lix_mat = zeros(Int8, NAx, 2)
    lix_border = zeros(Int8, NAx, 2)
    ind_row = zeros(Int64, 2*NAx+2*NAy+2*NAz)
    ind_col = zeros(Int64, 2*NAx+2*NAy+2*NAz)
    vals_A = zeros(Float64, 2*NAx+2*NAy+2*NAz)
    bars_Lp_x = zeros(Float64, NAx, 6)

    num_ele = 0

    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for cont in range(1, stop=Nx - 1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont+1][cont2][cont3] == 1))

                        pos = mapAx[From_3D_to_1D(cont, cont2, cont3, Nx - 1, Ny)] + 1
                        
                        ind_row[num_ele+1] = pos
                        ind_col[num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]], nodi_red)

                        bars_Lp_x[pos,1] = min_v[1] + sx * (cont - 1) + sx/2
                        bars_Lp_x[pos,2] = min_v[2] + sy * (cont2 - 1)
                        bars_Lp_x[pos,3] = min_v[3] + sz * (cont3 - 1)

                        vals_A[num_ele+1] = -1.0
                        num_ele = num_ele + 1

                        ind_row[num_ele+1] = pos
                        ind_col[num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont+1, cont2, cont3, Nx, Ny)]], nodi_red)
                        vals_A[num_ele+1] = 1.0
                        num_ele = num_ele + 1

                        bars_Lp_x[pos,4] = bars_Lp_x[pos,1] + sx
                        bars_Lp_x[pos,5] = bars_Lp_x[pos,2] + sy
                        bars_Lp_x[pos,6] = bars_Lp_x[pos,3] + sz

                        lix_mat[pos, 1] = k+1
                        lix_mat[pos, 2] = k+1

                        if cont > 1
                           if (matrice[k][cont-1][cont2][cont3]==0)
                               lix_border[pos, 1] = k + 1
                               bars_Lp_x[pos,1] = bars_Lp_x[pos,1] - sx/2.0
                           end
                        else
                            lix_border[pos, 1] = k + 1
                            bars_Lp_x[pos,1] = bars_Lp_x[pos,1] - sx / 2.0
                        end

                        if cont + 1 == Nx
                            lix_border[pos, 2] = k + 1
                            bars_Lp_x[pos,4] = bars_Lp_x[pos,4] + sx/2.0
                        else
                            if cont + 2 <= Nx
                                if (matrice[k][cont+2][cont2][cont3]==0)
                                    lix_border[pos, 2] = k + 1
                                    bars_Lp_x[pos,4] = bars_Lp_x[pos,4] + sx/2.0
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
    end


    liy_mat = zeros(Int8, NAy, 2)
    liy_border = zeros(Int8, NAy, 2)
    bars_Lp_y = zeros(Float64, NAy, 6)

    starter=num_ele
    num_ele = 0

    for cont3 in range(1, stop=Nz)
        for cont in range(1, stop=Nx)
            for cont2 in range(1, stop=Ny - 1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2+1][cont3] == 1))

                        pos = mapAy[From_3D_to_1D(cont, cont2, cont3, Nx, Ny - 1)] + 1

                        bars_Lp_y[pos,1] = min_v[1] + sx * (cont - 1)
                        bars_Lp_y[pos,2] = min_v[2] + sy * (cont2 - 1) + sy/2
                        bars_Lp_y[pos,3] = min_v[3] + sz * (cont3 - 1)

                        ind_row[starter+num_ele+1] = pos+NAx
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = -1.0
                        num_ele = num_ele + 1

                        bars_Lp_y[pos,4] = bars_Lp_y[pos,1] + sx
                        bars_Lp_y[pos,5] = bars_Lp_y[pos,2] + sy
                        bars_Lp_y[pos,6] = bars_Lp_y[pos,3] + sz

                        ind_row[starter+num_ele+1] = pos+NAx
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2+1, cont3, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = 1.0
                        num_ele = num_ele + 1

                        liy_mat[pos, 1] = k+1
                        liy_mat[pos, 2] = k+1

                        if cont2 > 1
                           if (matrice[k][cont][cont2-1][cont3]==0)
                               liy_border[pos, 1] = k + 1
                               bars_Lp_y[pos,2] = bars_Lp_y[pos,2] - sy/2.0
                           end
                        else
                            liy_border[pos, 1] = k + 1
                            bars_Lp_y[pos,2] = bars_Lp_y[pos,2] - sy / 2.0
                        end

                        if cont2 + 1 == Ny
                            liy_border[pos, 2] = k + 1
                            bars_Lp_y[pos,5] = bars_Lp_y[pos,5] + sy / 2.0
                        else
                            if cont2 + 2 <= Ny
                                if (matrice[k][cont][cont2+2][cont3]==0)
                                    liy_border[pos, 2] = k + 1
                                    bars_Lp_y[pos,5] = bars_Lp_y[pos,5] + sy / 2.0
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
    end


    liz_mat = zeros(Int8, NAz, 2)
    liz_border = zeros(Int8, NAz, 2)
    bars_Lp_z = zeros(Float64, NAz, 6)

    starter = starter+num_ele
    num_ele = 0

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz - 1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2][cont3+1] == 1))

                        pos = mapAz[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)] + 1

                        bars_Lp_z[pos,1] = min_v[1] + sx * (cont - 1)
                        bars_Lp_z[pos,2] = min_v[2] + sy * (cont2 - 1)
                        bars_Lp_z[pos,3] = min_v[3] + sz * (cont3 - 1) + sz/2

                        ind_row[starter+num_ele+1] = pos+NAx+NAy
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = -1.0
                        num_ele = num_ele + 1

                        bars_Lp_z[pos,4] = bars_Lp_z[pos,1] + sx
                        bars_Lp_z[pos,5] = bars_Lp_z[pos,2] + sy
                        bars_Lp_z[pos,6] = bars_Lp_z[pos,3] + sz

                        ind_row[starter+num_ele+1] = pos+NAx+NAy
                        ind_col[starter+num_ele+1] = bin_search(nodi[mapping_Vox[From_3D_to_1D(cont, cont2, cont3+1, Nx, Ny)]], nodi_red)
                        vals_A[starter+num_ele+1] = 1.0
                        num_ele = num_ele + 1

                        liz_mat[pos, 1] = k+1
                        liz_mat[pos, 2] = k+1

                        if cont3 > 1
                           if (matrice[k][cont][cont2][cont3-1]==0)
                               liz_border[pos, 1] = k + 1
                               bars_Lp_z[pos,3] = bars_Lp_z[pos,3] - sz/2.0
                           end
                        else
                            liz_border[pos, 1] = k + 1
                            bars_Lp_z[pos,3] = bars_Lp_z[pos,3] - sz / 2.0
                        end

                        if cont3 + 1 == Nz
                            liz_border[pos, 2] = k + 1
                            bars_Lp_z[pos,6] = bars_Lp_z[pos,6] + sz / 2.0
                        else
                            if cont3 + 2 <= Nz
                                if (matrice[k][cont][cont2][cont3+2]==0)
                                    liz_border[pos, 2] = k + 1
                                    bars_Lp_z[pos,6] = bars_Lp_z[pos,6] + sz / 2.0
                                end
                            end
                        end
                        break
                    end
                end
            end
        end
    end


    return ind_row,ind_col,vals_A,lix_mat,liy_mat,liz_mat,lix_border,liy_border,liz_border,bars_Lp_x,bars_Lp_y,bars_Lp_z
end

function ver_con(A,B)
    return vcat(A,B)
end

function vect_con(A,B)
    if ((typeof(A) == Vector{Int} && typeof(B) == Vector{Int}) || (typeof(A) == Vector{Float64} && typeof(B) == Vector{Float64}))
        return [A;B]
    else
        return hcat(A,B)
    end
end

function compute_diagonals(MATER,sx,sy,sz,lix_mat,liy_mat,liz_mat,lix_border,liy_border,liz_border)
    eps0 = 8.854187816997944e-12
    num_grids = length(MATER)
    for cont in range(1, stop=num_grids)
        sigmar = MATER[cont].conductivity
        epsr = MATER[cont].permittivity 
        if sigmar!=0
            MATER[cont].Rx = 0.5 * sx / (sigmar * sy * sz)
            MATER[cont].Ry = 0.5 * sy / (sigmar * sx * sz)
            MATER[cont].Rz = 0.5 * sz / (sigmar * sy * sx)
            if epsr == 1
                MATER[cont].Cx = 0.0
                MATER[cont].Cy = 0.0
                MATER[cont].Cz = 0.0
            else
                MATER[cont].Cx = eps0 * (epsr - 1.0) * sy * sz / (0.5 * sx)
                MATER[cont].Cy = eps0 * (epsr - 1.0) * sx * sz / (0.5 * sy)
                MATER[cont].Cz = eps0 * (epsr - 1.0) * sy * sx / (0.5 * sz)
            end
        else
            MATER[cont].Rx = 0.0
            MATER[cont].Ry = 0.0
            MATER[cont].Rz = 0.0
            MATER[cont].Cx = eps0 * (epsr - 1.0) * sy * sz / (0.5 * sx)
            MATER[cont].Cy = eps0 * (epsr - 1.0) * sx * sz / (0.5 * sy)
            MATER[cont].Cz = eps0 * (epsr - 1.0) * sy * sx / (0.5 * sz)
        end
    end
    Rx = zeros(Complex ,(size(lix_border)[1], 4))
    Ry = zeros(Complex ,(size(liy_border)[1], 4))
    Rz = zeros(Complex ,(size(liz_border)[1], 4))

    Cx = zeros(Complex ,(size(lix_border)[1], 4))
    Cy = zeros(Complex ,(size(liy_border)[1], 4))
    Cz = zeros(Complex ,(size(liz_border)[1], 4))

    for cont in range(1, stop=num_grids)
        if MATER[cont].Rx!=0
            ind_m = findall(l -> (cont+1) == l, lix_mat[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Rx[ind_m, 1] .= MATER[cont].Rx
            ind_m = findall(l -> (cont+1) == l, lix_mat[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Rx[ind_m, 2] .= MATER[cont].Rx
            ind_m = findall(l -> (cont+1) == l, lix_border[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Rx[ind_m, 3] .= MATER[cont].Rx
            ind_m = findall(l -> (cont+1) == l, lix_border[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Rx[ind_m, 4] .= MATER[cont].Rx
        end
        if MATER[cont].Cx!=0
            ind_m = findall(l -> (cont+1) == l, lix_mat[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Cx[ind_m, 1] .= MATER[cont].Cx
            ind_m = findall(l -> (cont+1) == l, lix_mat[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Cx[ind_m, 2] .= MATER[cont].Cx
            ind_m = findall(l -> (cont+1) == l, lix_border[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Cx[ind_m, 3] .= MATER[cont].Cx
            ind_m = findall(l -> (cont+1) == l, lix_border[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Cx[ind_m, 4] .= MATER[cont].Cx
        end
        if MATER[cont].Ry!=0
            ind_m = findall(l -> (cont+1) == l, liy_mat[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Ry[ind_m, 1] .= MATER[cont].Ry
            ind_m = findall(l -> (cont+1) == l, liy_mat[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Ry[ind_m, 2] .= MATER[cont].Ry
            ind_m = findall(l -> (cont+1) == l, liy_border[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Ry[ind_m, 3] .= MATER[cont].Ry
            ind_m = findall(l -> (cont+1) == l, liy_border[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Ry[ind_m, 4] .= MATER[cont].Ry
        end
        if MATER[cont].Cy!=0
            ind_m = findall(l -> (cont+1) == l, liy_mat[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Cy[ind_m, 1] .= MATER[cont].Cy
            ind_m = findall(l -> (cont+1) == l, liy_mat[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Cy[ind_m, 2] .= MATER[cont].Cy
            ind_m = findall(l -> (cont+1) == l, liy_border[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Cy[ind_m, 3] .= MATER[cont].Cy
            ind_m = findall(l -> (cont+1) == l, liy_border[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Cy[ind_m, 4] .= MATER[cont].Cy
        end
        if MATER[cont].Rz!=0
            ind_m = findall(l -> (cont+1) == l, liz_mat[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Rz[ind_m, 1] .= MATER[cont].Rz
            ind_m = findall(l -> (cont+1) == l, liz_mat[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Rz[ind_m, 2] .= MATER[cont].Rz
            ind_m = findall(l -> (cont+1) == l, liz_border[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Rz[ind_m, 3] .= MATER[cont].Rz
            ind_m = findall(l -> (cont+1) == l, liz_border[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Rz[ind_m, 4] .= MATER[cont].Rz
        end
        if MATER[cont].Cz!=0
            ind_m = findall(l -> (cont+1) == l, liz_mat[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Cz[ind_m, 1] .= MATER[cont].Cz
            ind_m = findall(l -> (cont+1) == l, liz_mat[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Cz[ind_m, 2] .= MATER[cont].Cz
            ind_m = findall(l -> (cont+1) == l, liz_border[:, 1])
            #ind_m = filter(i -> i!=1, ind_m)
            Cz[ind_m, 3] .= MATER[cont].Cz
            ind_m = findall(l -> (cont+1) == l, liz_border[:, 2])
            #ind_m = filter(i -> i!=1, ind_m)
            Cz[ind_m, 4] .= MATER[cont].Cz
        end
    end
    diag_R = ver_con(ver_con(Rx, Ry), Rz)
    diag_Cd = ver_con(ver_con(Cx, Cy), Cz)


    return diag_R, diag_Cd
end

function create_Gamma_and_center_sup(matrice, Nx,Ny,Nz, map_volumes, min_v, sx, sy, sz, nodi, nodi_red, ext_grids)
    num_grids = size(matrice)[1]
    mapping_surf_1 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_2 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_3 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_4 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_5 = zeros(Int64, Nx * Ny * Nz)
    mapping_surf_6 = zeros(Int64, Nx * Ny * Nz)

    num_ele_1 = 0
    num_ele_2 = 0
    num_ele_3 = 0
    num_ele_4 = 0
    num_ele_5 = 0
    num_ele_6 = 0

    nnz_surf_max = 6 * Nx * Ny * Nz
    ind_r = ones(Int64, nnz_surf_max)
    ind_c = ones(Int64, nnz_surf_max)

    sup_centers = zeros(Float64, nnz_surf_max,3)
    sup_type = zeros(Int64, nnz_surf_max)


    contat_tot = 1


    cont2=1
    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3]==1) && (ext_grids[k,2,cont,cont2,cont3]==0))
                    num_ele_1 = num_ele_1 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_1[p31] = num_ele_1
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1)  + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1)
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 1
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(2, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont][cont2-1][cont3]==0) && (ext_grids[k,2,cont,cont2,cont3]==0))
                        num_ele_1 = num_ele_1 + 1
                        p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        mapping_surf_1[p31] = num_ele_1
                        contat_tot = contat_tot + 1
                        ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                        ind_c[contat_tot] = num_ele_1
                        sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                        sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1)
                        sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                        sup_type[ind_c[contat_tot]] = 1
                        break
                    end
                end
            end
        end
    end

    cont2 = Ny
    for cont in range(1, stop=Nx)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,1,cont,cont2,cont3] == 0))
                    num_ele_2 = num_ele_2 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_2[p31] = num_ele_2
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 1
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny-1)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2+1][cont3] == 0) && (ext_grids[k,1,cont,cont2,cont3] == 0))
                        check_others = false
                        for k2 in range(1, stop=num_grids)
                            if k!=k2
                                if (matrice[k2][cont][cont2+1][cont3]==1)
                                    check_others = true
                                    break
                                end
                            end
                        end

                        if (check_others == false)
                            num_ele_2 = num_ele_2 + 1
                            p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            mapping_surf_2[p31] = num_ele_2
                            contat_tot = contat_tot + 1
                            ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                            ind_c[contat_tot] = num_ele_1 + num_ele_2
                            sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                            sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy
                            sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                            sup_type[ind_c[contat_tot]] = 1
                        end
                        break
                    end
                end
            end
        end
    end


    cont = 1
    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,4,cont,cont2,cont3] == 0))
                    num_ele_3 = num_ele_3 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_3[p31] = num_ele_3
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1)
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 2
                    break
                end
            end
        end
    end

    for cont in range(2, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont-1][cont2][cont3]==0) && (ext_grids[k,4,cont,cont2,cont3]==0))
                        num_ele_3 = num_ele_3 + 1
                        p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        mapping_surf_3[p31] = num_ele_3
                        contat_tot = contat_tot + 1
                        ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                        ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3
                        sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1)
                        sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                        sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                        sup_type[ind_c[contat_tot]] = 2
                        break
                    end
                end
            end
        end
    end

    cont = Nx
    for cont2 in range(1, stop=Ny)
        for cont3 in range(1, stop=Nz)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,3,cont,cont2,cont3] == 0))
                    num_ele_4 = num_ele_4 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_4[p31] = num_ele_4
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                    sup_type[ind_c[contat_tot]] = 2
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx-1)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont+1][cont2][cont3] == 0) && (ext_grids[k,3,cont,cont2,cont3] == 0))
                        check_others = false
                        for k2 in range(1, stop=num_grids)
                            if k!=k2
                                if (matrice[k2][cont+1][cont2][cont3]==1)
                                    check_others = true
                                    break
                                end
                            end
                        end

                        if (check_others == false)
                            num_ele_4 = num_ele_4 + 1
                            p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            mapping_surf_4[p31] = num_ele_4
                            contat_tot = contat_tot + 1
                            ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                            ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4
                            sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx
                            sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                            sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz / 2.0
                            sup_type[ind_c[contat_tot]] = 2
                        end
                        break
                    end
                end
            end
        end
    end



    cont3 = 1
    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,6,cont,cont2,cont3] == 0))
                    num_ele_5 = num_ele_5 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_5[p31] = num_ele_5
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1)
                    sup_type[ind_c[contat_tot]] = 3
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(2, stop=Nz)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3]==1) && (matrice[k][cont][cont2][cont3-1]==0) && (ext_grids[k,6,cont,cont2,cont3]==0))
                        num_ele_5 = num_ele_5 + 1
                        p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                        mapping_surf_5[p31] = num_ele_5
                        contat_tot = contat_tot + 1
                        ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                        ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5
                        sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                        sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                        sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1)
                        sup_type[ind_c[contat_tot]] = 3
                        break
                    end
                end
            end
        end
    end

    cont3 = Nz
    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for k in range(1, stop=num_grids)
                if ((matrice[k][cont][cont2][cont3] == 1) && (ext_grids[k,5,cont,cont2,cont3] == 0))
                    num_ele_6 = num_ele_6 + 1
                    p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                    mapping_surf_6[p31] = num_ele_6
                    contat_tot = contat_tot + 1
                    ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                    ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5 + num_ele_6
                    sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                    sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                    sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz
                    sup_type[ind_c[contat_tot]] = 3
                    break
                end
            end
        end
    end

    for cont in range(1, stop=Nx)
        for cont2 in range(1, stop=Ny)
            for cont3 in range(1, stop=Nz-1)
                for k in range(1, stop=num_grids)
                    if ((matrice[k][cont][cont2][cont3] == 1) && (matrice[k][cont][cont2][cont3+1] == 0) && (ext_grids[k,5,cont,cont2,cont3] == 0))
                        check_others = false
                        for k2 in range(1, stop=num_grids)
                            if k!=k2
                                if (matrice[k2][cont][cont2][cont3+1]==1)
                                    check_others = true
                                    break
                                end
                            end
                        end

                        if (check_others==false)
                            num_ele_6 = num_ele_6 + 1
                            p31 = From_3D_to_1D(cont, cont2, cont3, Nx, Ny)
                            mapping_surf_6[p31] = num_ele_6
                            contat_tot = contat_tot + 1
                            ind_r[contat_tot] = bin_search(nodi[map_volumes[p31]], nodi_red)
                            ind_c[contat_tot] = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5 + num_ele_6
                            sup_centers[ind_c[contat_tot-1], 1] = min_v[1] + sx * (cont - 1) + sx / 2.0
                            sup_centers[ind_c[contat_tot-1], 2] = min_v[2] + sy * (cont2 - 1) + sy / 2.0
                            sup_centers[ind_c[contat_tot-1], 3] = min_v[3] + sz * (cont3 - 1) + sz
                            sup_type[ind_c[contat_tot]] = 3
                        end
                        break
                    end
                end
            end
        end
    end

    nnz_surf = num_ele_1 + num_ele_2 + num_ele_3 + num_ele_4 + num_ele_5 + num_ele_6
    vals = ones(Float64, contat_tot)
    ind_r = ind_r[1:contat_tot]
    ind_c = ind_c[1:contat_tot]
    sup_centers = sup_centers[1:nnz_surf,:]
    sup_type = sup_type[1:nnz_surf]


    return vals,ind_r,ind_c,sup_centers,sup_type
end

function generate_interconnection_matrices_and_centers(size_x,size_y,size_z,grid_matrix,num_cel_x,num_cel_y,num_cel_z,materials,port_matrix,lumped_el_matrix,minimum_vertex)

    numTotVox = num_cel_x*num_cel_y*num_cel_z
    
    println("Total Number of Voxels (including air):", numTotVox)
    n_grids = length(materials)
    @assert size(grid_matrix)[1]==n_grids
    
    num_tot_full_vox = 0
    
    for i in range(1, stop=n_grids)
        for j in range(1, stop=num_cel_x)
            for k in range(1, stop=num_cel_y)
                num_tot_full_vox = num_tot_full_vox + count(i-> i==1,grid_matrix[i][j][k])
            end
        end
    end

    # TODO: VERIFICARE SE OCCORRE
    cont=1
    while cont<=length(materials)
        materials[cont].epsr = materials[cont].permittivity +1im * materials[cont].permittivity * materials[cont].tangent_delta_permittivity
        cont+=1
    end

    @assert num_cel_x isa Int
    @assert num_cel_y isa Int
    @assert num_cel_z isa Int

    n_boxes,mapping_vols,volume_centers,volumes_materials=create_volumes_mapping_and_centers(grid_matrix,num_cel_x,num_cel_y,num_cel_z,num_tot_full_vox,size_x,size_y,size_z,minimum_vertex)

    externals_grids = create_external_grids(grid_matrix,num_cel_x,num_cel_y,num_cel_z)
    
    nodes_red, nodes = create_nodes_ref(grid_matrix,num_cel_x,num_cel_y,num_cel_z,num_tot_full_vox,externals_grids,mapping_vols)

    port_matrix.port_voxels, port_matrix.port_nodes = find_nodes_port(volume_centers, port_matrix.port_start, port_matrix.port_end, nodes, nodes_red)
        

    lumped_el_matrix.le_voxels, lumped_el_matrix.le_nodes = find_nodes_port(volume_centers,lumped_el_matrix.le_start,lumped_el_matrix.le_end, nodes, nodes_red)

    vals,ind_r,ind_c,sup_centers,sup_type=create_Gamma_and_center_sup(grid_matrix,num_cel_x,num_cel_y,num_cel_z,mapping_vols, minimum_vertex, size_x,size_y,size_z, nodes, nodes_red, externals_grids)

    Gamma = sparse(ind_r, ind_c, vals)

    println("Number of Surfaces (without air):", size(Gamma)[2])

    map_for_Ax, n_for_Ax = create_mapping_Ax(grid_matrix,num_cel_x,num_cel_y,num_cel_z)
    map_for_Ay, n_for_Ay = create_mapping_Ay(grid_matrix,num_cel_x,num_cel_y,num_cel_z)
    map_for_Az, n_for_Az = create_mapping_Az(grid_matrix,num_cel_x,num_cel_y,num_cel_z)

    ind_row,ind_col,vals_A, LiX, LiY, LiZ, LiX_bord,LiY_bord,LiZ_bord,bars_Lp_x,bars_Lp_y,bars_Lp_z = create_A_mats_volInd(grid_matrix,num_cel_x,num_cel_y,num_cel_z,mapping_vols,
                                       map_for_Ax, n_for_Ax, map_for_Ay, n_for_Ay, map_for_Az, n_for_Az,
                                       size_x,size_y,size_z,minimum_vertex,nodes,nodes_red)

    

    A = sparse(ind_row, ind_col, vals_A)

    println("Edges without air:", n_for_Ax+n_for_Ay+n_for_Az)
    println("Nodes without air:", size(Gamma)[1])

    diag_R,diag_Cd=compute_diagonals(materials, size_x,size_y,size_z, LiX, LiY, LiZ,LiX_bord, LiY_bord, LiZ_bord)
    

    return A, Gamma, port_matrix, lumped_el_matrix, sup_centers, sup_type, bars_Lp_x, bars_Lp_y, bars_Lp_z, diag_R, diag_Cd
end