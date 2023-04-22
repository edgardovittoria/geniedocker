using FLoops
function compute_Lp_self(l,W,T)
    #fast Henry
    w=W/l
    t=T/l
    r = sqrt(w * w + t * t)
    aw = sqrt(w * w + 1.0)
    at = sqrt(t * t + 1.0)
    ar = sqrt(w * w + t * t + 1.0)

    mu0 = 4.0 * pi * 1e-7

    Lp_Self_Rect = 2.0 * mu0 * l / pi * ( 0.25 * (1 / w * asinh(w / at) + 
                    1 / t * asinh(t / aw) + asinh(1.0 / r)) + 
                1.0 / 24.0 * (t * t / w * asinh(w / (t * at * (r + ar))) 
                              + w * w / t * asinh(t / (w * aw * (r + ar))) + 
                              t * t / (w * w) * asinh(w * w / (t * r * (at + ar))) + 
                              w * w / (t * t) * asinh(
                    t *t / (w * r * (aw + ar))) + 
                              1.0 / (w * t * t) * asinh(w * t * t / (at * (aw + ar))) + 1.0 / (
                                          t * w * w) * asinh(t * w * w / (aw * (at + ar))) 
                              ) 
                - 1.0 / 6.0 * (
                            1.0 / (w * t) * atan(w * t / ar) + t / w * atan(w / (t * ar)) 
                            + w / t * atan(t / (w * ar))) 
                - 1.0 / 60.0 * ((ar + r + t + at) * t * t / ((ar + r) * (r + t) * (t + at) * (at + ar)) + 
                                (ar + r + w + aw) * w * w / ((ar + r) * (r + w) * (w + aw) * (aw + ar)) + 
                                (ar + aw + 1.0 + at) / ((ar + aw) * (aw + 1.0) * (1.0 + at) * (at + ar)) 
                                ) 
                - 1.0 / 20.0 * (1.0 / (r + ar) + 1.0 / (aw + ar) + 1.0 / (at + ar))
                )
    return Lp_Self_Rect
end


function Integ_vol_vol(x1v,y1v,z1v, x2v,y2v,z2v)

    epsilon1=5e-3
    epsilon2=1e-3
    epsilon3=1e-3
    epsilon4=3e-1

    xc1=0.5*(x1v[1]+x1v[2])
    yc1=0.5*(y1v[1]+y1v[2])
    zc1=0.5*(z1v[1]+z1v[2])
    xc2=0.5*(x2v[1]+x2v[2])
    yc2=0.5*(y2v[1]+y2v[2])
    zc2=0.5*(z2v[1]+z2v[2])

    a1=abs(x1v[2]-x1v[1])
    b1=abs(y1v[2]-y1v[1])
    c1=abs(z1v[2]-z1v[1])

    a2=abs(x2v[2]-x2v[1])
    b2=abs(y2v[2]-y2v[1])
    c2=abs(z2v[2]-z2v[1])

    V1=a1*b1*c1
    V2=a2*b2*c2

    supp_x1=false
    supp_y1=false
    supp_z1=false
    supp_x2=false
    supp_y2=false
    supp_z2=false

    aux_x = [abs(x1v[1] - x2v[1]), abs(x1v[1] - x2v[2]), abs(x1v[2] - x2v[1]), abs(x1v[2] - x2v[2])]
    aux_y = [abs(y1v[1] - y2v[1]), abs(y1v[1] - y2v[2]), abs(y1v[2] - y2v[1]), abs(y1v[2] - y2v[2])]
    aux_z = [abs(z1v[1] - z2v[1]), abs(z1v[1] - z2v[2]), abs(z1v[2] - z2v[1]), abs(z1v[2] - z2v[2])]

    min_R = sqrt((minimum(aux_x)^2) + (minimum(aux_y)^2) + (minimum(aux_z)^2))

    if (a1 <= b1 && a1 <= c1)
        max_d = maximum(aux_x)
        supp_x1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, a1, b1, c1)
        if supp_x1 == true
            max_ed = max(b1, c1)
            if (max_ed / (min_R + 1e-15) < epsilon4)
                max_d = maximum(aux_y)
                supp_y1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, b1, a1, c1)
                max_d = maximum(aux_z)
                supp_z1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, c1, a1, b1)
            end
        end
    else
        if (b1 <= a1 && b1 <= c1)
            max_d = maximum(aux_y)
            supp_y1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, b1, a1, c1)
            if supp_y1 == true
                max_ed = max(a1, c1)
                if (max_ed / (min_R + 1e-15) < epsilon4)
                    max_d = maximum(aux_x)
                    supp_x1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, a1, b1, c1)
                    max_d = maximum(aux_z)
                    supp_z1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, c1, a1, b1)
                end
            end
        else
            max_d = maximum(aux_z)
            supp_z1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, c1, a1, b1)
            if supp_z1 == true
                max_ed = max(a1, b1)
                if (max_ed / (min_R + 1e-15) < epsilon4)
                    max_d = maximum(aux_x)
                    supp_x1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, a1, b1, c1)
                    max_d = maximum(aux_y)
                    supp_y1 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, b1, a1, c1)
                end
            end
        end
    end

    if (a2 <= b2 && a2 <= c2)
        max_d = maximum(aux_x)
        supp_x2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, a2, b2, c2)
        if supp_x2 == true
            max_ed = max(b2, c2)
            if (max_ed / (min_R + 1e-15) < epsilon4)
                max_d = maximum(aux_y)
                supp_y2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, b2, a2, c2)
                max_d = maximum(aux_z)
                supp_z2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, c2, a2, b2)
            end
        end
    else
        if (b2 <= a2 && b2 <= c2)
            max_d = maximum(aux_y)
            supp_y2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, b2, a2, c2)
            if supp_y2 == true
                max_ed = max(a2, c2)
                if (max_ed / (min_R + 1e-15) < epsilon4)
                    max_d = maximum(aux_x)
                    supp_x2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, a2, b2, c2)
                    max_d = maximum(aux_z)
                    supp_z2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, c2, a2, b2)
                end
            end
        else
            max_d = maximum(aux_z)
            supp_z2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, c2, a2, b2)
            if supp_z2 == true
                max_ed = max(a2, b2)
                if (max_ed / (min_R + 1e-15) < epsilon4)
                    max_d = maximum(aux_x)
                    supp_x2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, a2, b2, c2)
                    max_d = maximum(aux_y)
                    supp_y2 = check_condition(epsilon1, epsilon2, epsilon3, V1, V2, max_d, min_R, b2, a2, c2)
                end
            end
        end
    end

    sum_supp = supp_x1 + supp_y1 + supp_z1 + supp_x2 + supp_y2 + supp_z2

    # point - point to use
    if sum_supp == 6
        integ = V1 * V2 / sqrt( ((xc1 - xc2)^2) + ((yc1 - yc2)^2) + ((zc1 - zc2)^2) )
    elseif sum_supp == 5  #point-line to use
        is_point_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 3)
            is_point_v1 = true
        end

        if is_point_v1 == true
            if supp_x2 == false # line of volume 2 along x
                integ = V1 * V2 / a2 * integ_line_point(x2v, yc2, zc2, xc1, yc1, zc1)
            else
                if supp_y2 == false #line of volume 2 along y
                    integ = V1 * V2 / b2 * integ_line_point(y2v, xc2, zc2, yc1, xc1, zc1)
                else # line of volume 2 along z
                    integ = V1 * V2 / c2 * integ_line_point(z2v, xc2, yc2, zc1, xc1, yc1)
                end
            end
        else
            if supp_x1 == false # line of volume 1 along x
                integ = V1 * V2 / a1 * integ_line_point(x1v, yc1, zc1, xc2, yc2, zc2)
            else
                if supp_y1 == false # line of volume 1 along y
                    integ = V1 * V2 / b1 * integ_line_point(y1v, xc1, zc1, yc2, xc2, zc2)
                else # %line of volume 1 along z
                    integ = V1 * V2 / c1 * integ_line_point(z1v, xc1, yc1, zc2, xc2, yc2)
                end
            end
        end
    elseif sum_supp == 4 #point-surface or line-line case
        is_point_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 3)
            is_point_v1 = true
        end
        is_point_v2 = false
        if (supp_x2 + supp_y2 + supp_z2 == 3)
            is_point_v2 = true
        end

        if is_point_v1 == true # point - surface case
            if supp_x2 == true # surface of volume 2 in yz plane
                integ = V1 * a2 * integ_point_sup(zc1, yc1, xc1, z2v, y2v, xc2)
            else
                if supp_y2 == true #surface of volume 2 in xz plane
                    integ = V1 * b2 * integ_point_sup(xc1, zc1, yc1, x2v, z2v, yc2)
                else #surface of volume 2 in xy plane
                    integ = V1 * c2 * integ_point_sup(xc1, yc1, zc1, x2v, y2v, zc2)
                end
            end
        else
            if is_point_v2 == true #point-surface case
                if supp_x1 == true #surface of volume 1 in yz plane
                    integ = V2 * a1 * integ_point_sup(zc2, yc2, xc2, z1v, y1v, xc1)
                else
                    if supp_y1 == true #surface of volume 1 in xz plane
                        integ = V2 * b1 * integ_point_sup(xc2, zc2, yc2, x1v, z1v, yc1)
                    else #surface of volume 1 in xy plane
                        integ = V2 * c1 * integ_point_sup(xc2, yc2, zc2, x1v, y1v, zc1)
                    end
                end
            else #line-line case
                if supp_y1 == true && supp_z1 == true
                    if supp_y2 == true && supp_z2 == true # parallel lines
                        integ = b1 * c1 * b2 * c2 * integ_line_line_parall(x1v, yc1, zc1, x2v, yc2, zc2)
                    else
                        if supp_x2 == true && supp_z2 == true # orthogonal lines
                            integ = b1 * c1 * a2 * c2 * integ_line_line_ortho_xy(x1v, yc1, zc1, xc2, y2v, zc2)
                        else
                            integ = b1 * c1 * a2 * b2 * integ_line_line_ortho_xy(x1v, zc1, yc1, xc2, z2v, yc2)
                        end
                    end
                else
                    if supp_x1 == true && supp_z1 == true
                        if supp_x2 == true && supp_z2 == true #parallel lines
                            integ = a1 * c1 * a2 * c2 * integ_line_line_parall(y1v, xc1, zc1, y2v, xc2, zc2)
                        else
                            if supp_x2 == true && supp_y2 == true # orthogonal lines
                                integ = a1 * c1 * a2 * b2 * integ_line_line_ortho_xy(y1v, zc1, xc1, yc2, z2v, xc2)
                            else
                                integ = a1 * c1 * b2 * c2 * integ_line_line_ortho_xy(y1v, xc1, zc1, yc2, x2v, zc2)
                            end
                        end
                    else
                        if supp_x2 == true && supp_y2 == true # parallel lines
                            integ = a1 * b1 * a2 * b2 * integ_line_line_parall(z1v, xc1, yc1, z2v, xc2, yc2)
                        else
                            if supp_x2 == true && supp_z2 == true # orthogonal lines
                                integ = a1 * b1 * a2 * c2 * integ_line_line_ortho_xy(z1v, yc1, xc1, zc2, y2v, xc2)
                            else
                                integ = a1 * b1 * b2 * c2 * integ_line_line_ortho_xy(z1v, xc1, yc1, zc2, x2v, yc2)
                            end
                        end
                    end
                end
            end
        end
    elseif sum_supp == 3  # point-volume or surface-line
        is_point_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 3)
            is_point_v1 = true
        end
        is_point_v2 = false
        if (supp_x2 + supp_y2 + supp_z2 == 3)
            is_point_v2 = true
        end
        is_surf_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 1)
            is_surf_v1 = true
        end
        if is_point_v1 == true # point - volume case
            integ = a1 * b1 * c1 * integ_point_vol(xc1, yc1, zc1, x2v, y2v, z2v)
        else
            if is_point_v2 == true # point - volume case
                integ = a2 * b2 * c2 * integ_point_vol(xc2, yc2, zc2, x1v, y1v, z1v)
            else # line-surface case
                if is_surf_v1 == true  # bar1 is a surface
                    if supp_x1 == true  # bar1 is a surface in y-z plane
                        if supp_x2 == false  # bar 2 is a line along x
                            integ = a1 * b2 * c2 * integ_line_surf_ortho(x2v, yc2, zc2, xc1, y1v, z1v)
                        else
                            if supp_y2 == false  # bar 2 is a line along y
                                integ = a1 * a2 * c2 * integ_line_surf_para(y1v, z1v, xc1, y2v, zc2, xc2)
                            else  # bar 2 is a line along z
                                integ = a1 * a2 * b2 * integ_line_surf_para(z1v, y1v, xc1, z2v, yc2, xc2)
                            end
                        end
                    else
                        if supp_y1 == true  # bar1 is a surface in x-z plane
                            if supp_x2 == false  # bar 2 is a line along x
                                integ = b1 * b2 * c2 * integ_line_surf_para(x1v, z1v, yc1, x2v, zc2, yc2)
                            else
                                if supp_y2 == false  # bar 2 is a line along y
                                    integ = b1 * a2 * c2 * integ_line_surf_ortho(y2v, xc2, zc2, yc1, x1v, z1v)
                                else  # bar 2 is a line along z
                                    integ = b1 * a2 * b2 * integ_line_surf_para(z1v, x1v, yc1, z2v, xc2, yc2)
                                end
                            end
                        else # bar1 is a surface in x-y plane
                            if supp_x2 == false  # bar 2 is a line along x
                                integ = c1 * b2 * c2 * integ_line_surf_para(x1v, y1v, zc1, x2v, yc2, zc2)
                            else
                                if supp_y2 == false  # bar 2 is a line along y
                                    integ = c1 * a2 * c2 * integ_line_surf_para(y1v, x1v, zc1, y2v, xc2, zc2)
                                else  # bar 2 is a line along z
                                    integ = c1 * a2 * b2 * integ_line_surf_ortho(z2v, xc2, yc2, zc1, x1v, y1v)
                                end
                            end
                        end
                    end
                else  # bar2 is a surface
                    if supp_x2 == true  # bar2 is a surface in y-z plane
                        if supp_x1 == false  # bar 1 is a line along x
                            integ = a2 * b1 * c1 * integ_line_surf_ortho(x1v, yc1, zc1, xc2, y2v, z2v)
                        else
                            if supp_y1 == false  # bar 1 is a line along y
                                integ = a2 * a1 * c1 * integ_line_surf_para(y2v, z2v, xc2, y1v, zc1, xc1)
                            else  # bar 1 is a line along z
                                integ = a2 * a1 * b1 * integ_line_surf_para(z2v, y2v, xc2, z1v, yc1, xc1)
                            end
                        end
                    else
                        if supp_y2 == true  # bar2 is a surface in x-z plane
                            if supp_x1 == false  # bar 1 is a line along x
                                integ = b2 * b1 * c1 * integ_line_surf_para(x2v, z2v, yc2, x1v, zc1, yc1)
                            else
                                if supp_y1 == false  # bar 1 is a line along y
                                    integ = b2 * a1 * c1 * integ_line_surf_ortho(y1v, xc1, zc1, yc2, x2v, z2v)
                                else # bar 1 is a line along z
                                    integ = b2 * a1 * b1 * integ_line_surf_para(z2v, x2v, yc2, z1v, xc1, yc1)
                                end
                            end
                        else  # bar2 is a surface in x-y plane
                            if supp_x1 == false  # bar 1 is a line along x
                                integ = c2 * b1 * c1 * integ_line_surf_para(x2v, y2v, zc2, x1v, yc1, zc1)
                            else
                                if supp_y1 == false  # bar 1 is a line along y
                                    integ = c2 * a1 * c1 * integ_line_surf_para(y2v, x2v, zc2, y1v, xc1, zc1)
                                else  # bar 1 is a line along z
                                    integ = c2 * a1 * b1 * integ_line_surf_ortho(z1v, xc1, yc1, zc2, x2v, y2v)
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif sum_supp == 2 # line-volume or surface-surface
        is_line_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 2)
            is_line_v1 = true
        end
        is_line_v2 = false
        if (supp_x2 + supp_y2 + supp_z2 == 2)
            is_line_v2 = true
        end
        if is_line_v1 == true  # bar1 is a line
            if supp_x1 == false  # bar1 is a line along x
                integ = b1 * c1 * integ_line_vol(x2v, y2v, z2v, x1v, yc1, zc1)
            else
                if supp_y1 == false  # bar1 is a line along y
                    integ = a1 * c1 * integ_line_vol(y2v, x2v, z2v, y1v, xc1, zc1)
                else  # bar1 is a line along z
                    integ = a1 * b1 * integ_line_vol(z2v, x2v, y2v, z1v, xc1, yc1)
                end
            end
        else
            if is_line_v2 == true  # bar2 is a line
                if supp_x2 == false  # bar2 is a line along x
                    integ = b2 * c2 * integ_line_vol(x1v, y1v, z1v, x2v, yc2, zc2)
                else
                    if supp_y2 == false  # bar2 is a line along y
                        integ = a2 * c2 * integ_line_vol(y1v, x1v, z1v, y2v, xc2, zc2)
                    else  # bar2 is a line along z
                        integ = a2 * b2 * integ_line_vol(z1v, x1v, y1v, z2v, xc2, yc2)
                    end
                end
            else  # surface-surface case
                if supp_x1 == true  # bar1 is a surface in yz plane
                    if supp_x2 == true  # bar2 is a surface in yz plane
                        integ = a1 * a2 * integ_surf_surf_para(y1v, z1v, xc1, y2v, z2v, xc2)
                    else
                        if supp_y2 == true  # bar2 is a surface in xz plane
                            integ = a1 * b2 * integ_surf_surf_ortho(z1v, y1v, xc1, z2v, yc2, x2v)
                        else  # bar2 is a surface in xy plane
                            integ = a1 * c2 * integ_surf_surf_ortho(y1v, z1v, xc1, y2v, zc2, x2v)
                        end
                    end
                else
                    if supp_y1 == true  # bar1 is a surface in xz plane
                        if supp_x2 == true  # bar2 is a surface in yz plane
                            integ = b1 * a2 * integ_surf_surf_ortho(z1v, x1v, yc1, z2v, xc2, y2v)
                        else
                            if supp_y2 == true  # bar2 is a surface in xz plane
                                integ = b1 * b2 * integ_surf_surf_para(x1v, z1v, yc1, x2v, z2v, yc2)
                            else  # bar2 is a surface in xy plane
                                integ = b1 * c2 * integ_surf_surf_ortho(x1v, z1v, yc1, x2v, zc2, y2v)
                            end
                        end
                    else  # bar1 is a surface in xy plane

                        if supp_x2 == true  # bar2 is a surface in yz plane
                            integ = c1 * a2 * integ_surf_surf_ortho(y1v, x1v, zc1, y2v, xc2, z2v)
                        else
                            if supp_y2 == true  # bar2 is a surface in xz plane
                                integ = c1 * b2 * integ_surf_surf_ortho(x1v, y1v, zc1, x2v, yc2, z2v)
                            else  # bar2 is a surface in xy plane
                                integ = c1 * c2 * integ_surf_surf_para(x1v, y1v, zc1, x2v, y2v, zc2)
                            end
                        end
                    end
                end
            end
        end
    elseif sum_supp == 1  # surface-volume case
        if supp_x1 == true  # bar1 is a surface in yz plane
            integ = a1 * integ_vol_surf(y2v, z2v, x2v, y1v, z1v, xc1)
        else
            if supp_y1 == true  # bar1 is a surface in xz plane
                integ = b1 * integ_vol_surf(x2v, z2v, y2v, x1v, z1v, yc1)
            else
                if supp_z1 == true  # bar1 is a surface in xy plane
                    integ = c1 * integ_vol_surf(x2v, y2v, z2v, x1v, y1v, zc1)
                else
                    if supp_x2 == true  # bar2 is a surface in yz plane
                        integ = a2 * integ_vol_surf(y1v, z1v, x1v, y2v, z2v, xc2)
                    else
                        if supp_y2 == true  # bar2 is a surface in xz plane
                            integ = b2 * integ_vol_surf(x1v, z1v, y1v, x2v, z2v, yc2)
                        else # bar2 is a surface in xy plane
                            integ = c2 * integ_vol_surf(x1v, y1v, z1v, x2v, y2v, zc2)
                        end
                    end
                end
            end
        end
    else  # volume - volume case
        integ = integ_vol_vol(x1v, y1v, z1v, x2v, y2v, z2v)
    end
    return integ
end


function check_condition(eps1,eps2,eps3,V1,V2,max_d,min_R,size_dim,other_dim1,other_dim2)
    max_oth=maximum([other_dim1,other_dim2])
    condX1a = V1*V2*max_d/((min_R+1e-15)^3)
    condX1f = size_dim/(min_R+1e-15)
    condX1b = size_dim/max_oth
    supp_dim=false
    if ((condX1b <= eps3 || condX1f < eps1) && condX1a < eps2)
        supp_dim=true
    end
    return supp_dim
end


function integ_line_point(x1v,y3,z3,x1,y1,z1)
    x3=x1v[0]
    x4=x1v[1]

    if (x1 - x3)==0
        x1 = x1 - 1e-8
    end

    if (x1 - x4)==0
        x1 = x1 + 1e-8
    end

    R1 = sqrt(((x1 - x3)^2) + ((y1 - y3)^2) + ((z1 - z3)^2))
    R2 = sqrt(((x1 - x4)^2) + ((y1 - y3)^2) + ((z1 - z3)^2))

    Ip = real(log(complex((x1-x3)+R1+1e-20))) - real(log(complex((x1-x4)+R2+1e-20)))

    if (isnan(Ip) == true || isinf(Ip) == true)
        Ip = (x1 - x3) / (1e-14*abs((x1 - x3))) * real(log(complex(abs(x1-x3)+1e-20)))
             - (x1 - x4) / (1e-14*abs((x1 - x4))) * real(log(complex(abs(x1-x4)+1e-20)))
    end

    return Ip
end


function integ_point_sup(x1,y1,z1,x2v,y2v,z2)
    sol = 0
    for c1 in range(1, stop=2)
        x2 = x2v[c1]
        for c2 in range(1, stop=2)
            y2 = y2v[c2]

            R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

            if abs((y1 - y2) + R) < 1e-16
                term1 = 0.0
            else
                term1 = (x1 - x2) * real(log(complex((y1-y2)+R)))
            end

            if abs((x1 - x2) + R) < 1e-16
                term2 = 0.0
            else
                term2 = (y1 - y2) * real(log(complex((x1-x2)+R)))
            end

            term3 = -abs(z1 - z2) * atan((x1 - x2) * (y1 - y2), (abs(z1 - z2) * R))

            sol = sol + ((-1)^(c1 + c2 - 2)) * (term1 + term2 + term3)
        end
    end

    return sol
end


function integ_line_line_parall(x1v,y1,z1,x2v,y2,z2)
    sol = 0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            x2 = x2v[c2]

            R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

            if abs((x1 - x2) + R) < 1e-16
                term1 = 0.0
            else
                term1 = (x1 - x2) * real(log(complex((x1-x2)+R)))
            end

            sol = sol + ((-1)^(c1 + c2 + 1 - 2)) * (term1 - R)
        end
    end
    return sol
end


function integ_line_line_ortho_xy(x1v,y1,z1,x2,y2v,z2)
    sol = 0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y2 = y2v[c2]

            R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

            if abs((y1 - y2) + R) < 1e-16
                term1 = 0.0
            else
                term1 = (x1 - x2) * real(log(complex((y1-y2)+R)))
            end

            if abs((x1 - x2) + R) < 1e-16
                term2 = 0.0
            else
                term2 = (y1 - y2) * real(log(complex((x1-x2)+R)))
            end

            term3 = -abs(z1 - z2) * atan((x1 - x2) * (y1 - y2), (abs(z1 - z2) * R))

            sol = sol + ((-1)^(c1 + c2 + 1 - 2)) * (term1 + term2 + term3)
        end
    end

    return sol
end


function integ_point_vol(x1,y1,z1, x2v,y2v,z2v)
    sol = 0.0
    for c1 in range(1, stop=2)
        x2 = x2v[c1]
        for c2 in range(1, stop=2)
            y2 = y2v[c2]
            for c3 in range(1, stop=2)
                z2 = z2v[c3]
                R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                if abs((x1 - x2) + R) < 1e-16
                    term1 = 0.0
                else
                    term1 = (y1 - y2) * (z1 - z2) * real(log(complex((x1 - x2) + R)))
                end

                if abs((y1 - y2) + R) < 1e-16
                    term2 = 0.0
                else
                    term2 = (x1 - x2) * (z1 - z2) * real(log(complex((y1 - y2) + R)))
                end

                if abs((z1 - z2) + R) < 1e-16
                    term3 = 0.0
                else
                    term3 = (y1 - y2) * (x1 - x2) * real(log(complex((z1 - z2) + R)))
                end


                term4 = -0.5 * abs(z1 - z2) * (z1 - z2) * atan((x1 - x2) * (y1 - y2), (abs(z1 - z2) * R))

                term5 = -0.5 * abs(y1 - y2) * (y1 - y2) * atan((x1 - x2) * (z1 - z2), (abs(y1 - y2) * R))

                term6 = -0.5 * abs(x1 - x2) * (x1 - x2) * atan((y1 - y2) * (z1 - z2), (abs(x1 - x2) * R))

                sol = sol + ((-1.0)^(c1 + c2 + c3 - 3)) * (term1 + term2 + term3 + term4 + term5 + term6)
            end
        end
    end

    return sol
end


function integ_line_surf_ortho(x1v,y1,z1, x2,y2v,z2v)
    sol = 0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y2 = y2v[c2]
            for c3 in range(1, stop=2)
                z2 = z2v[c3]

                R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                if abs((x1 - x2) + R) < 1e-16
                    term1 = 0.0
                else
                    term1 = (y1 - y2) * (z1 - z2) * real(log(complex((x1-x2)+R)))
                end

                if abs((y1 - y2) + R) < 1e-16
                    term2 = 0.0
                else
                    term2 = (x1 - x2) * (z1 - z2) * real(log(complex((y1-y2))))
                end

                if abs((z1 - z2) + R) < 1e-16
                    term3 = 0.0
                else
                    term3 = (y1 - y2) * (x1 - x2) * real(log(complex((z1-z2)+R)))
                end

                term4 = -0.5 * abs(z1 - z2) * (z1 - z2) * atan((x1 - x2) * (y1 - y2), (abs(z1 - z2) * R))

                term5 = -0.5 * abs(y1 - y2) * (y1 - y2) * atan((x1 - x2) * (z1 - z2), (abs(y1 - y2) * R))

                term6 = -0.5 * abs(x1 - x2) * (x1 - x2) * atan((y1 - y2) * (z1 - z2), (abs(x1 - x2) * R))

                sol = sol + ((-1.0)^(c1 + c2 + c3+1 - 3))*(term1 + term2 + term3 + term4 + term5 + term6)
            end
        end
    end

    return sol
end


function integ_line_surf_para( x1v,y1v,z1,x2v,y2,z2)
    sol = 0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y1 = y1v[c2]
            for c3 in range(1, stop=2)
                x2 = x2v[c3]

                R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                if abs((x1 - x2) + R) < 1e-16
                    term1 = 0.0
                else
                    term1 = (x1 - x2) * (y1 - y2) * real(log(complex((x1-x2)+R)))
                end

                if abs((y1 - y2) + R) < 1e-16
                    term2 = 0.0
                else
                    term2 = (((x1 - x2)^2) - ((z1 - z2)^2)) / 2.0 * real(log(complex((y1-y2)+R)))
                end

                term3 = -(x1 - x2) * abs(z1 - z2) * atan((x1 - x2) * (y1 - y2), (abs(z1 - z2) * R))

                term4 = -(y1 - y2) / 2.0 * R

                sol = sol + ((-1.0)^(c1 + c2 + c3 - 3)) * (term1 + term2 + term3 + term4)
            end
        end
    end

    return sol
end


function integ_line_vol(x1v,y1v,z1v,x2v,y2,z2)
    sol = 0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y1 = y1v[c2]
            for c3 in range(1, stop=2)
                z1 = z1v[c3]
                for c4 in range(1, stop=2)
                    x2 = x2v[c4]
                    R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                    term1 = -1.0 / 3.0 * (y1 - y2) * (z1 - z2) * R

                    term2 = -0.5 * (x1 - x2) * abs(z1 - z2) * (z1 - z2) * atan((x1 - x2) * (y1 - y2),
                                                                                  (abs(z1 - z2) * R))

                    term3 = -0.5 * (x1 - x2) * abs(y1 - y2) * (y1 - y2) * atan((x1 - x2) * (z1 - z2),
                                                                                  (abs(y1 - y2) * R))

                    term4 = -1.0 / 6.0 * abs(((x1 - x2)^3)) * atan((y1 - y2) * (z1 - z2), (abs(x1 - x2) * R))

                    if abs((x1 - x2) + R) < 1e-16
                        term5 = 0.0
                    else
                        term5 = (x1 - x2) * (y1 - y2) * (z1 - z2) * real(log(complex((x1-x2)+R)))
                    end

                    if abs((y1 - y2) + R) < 1e-16
                        term6 = 0.0
                    else
                        term6 = (0.5 * ((x1 - x2)^2) - 1.0 / 6.0 * ((z1 - z2)^2)) * (z1 - z2) *
                            real(log(complex((y1-y2)+R)))
                    end

                    if abs((z1 - z2) + R) < 1e-16
                        term7 = 0.0
                    else
                        term7 = (0.5 * ((x1 - x2)^2) - 1.0 / 6.0 * ((y1 - y2)^2) ) * (y1 - y2) *
                            real(log(complex((z1 - z2) + R)))
                    end


                    sol = sol + ((-1.0)^(c1 + c2 + c3 + c4 + 1 - 4)) * (
                                term1 + term2 + term3 + term4 + term5 + term6 + term7)
                end
            end
        end
    end

    return sol
end


function integ_surf_surf_para(x1v,y1v,z1, x2v,y2v,z2)
    sol=0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y1 = y1v[c2]
            for c3 in range(1, stop=2)
                y2 = y2v[c3]
                for c4 in range(1, stop=2)
                    x2 = x2v[c4]
                    R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                    term1 = -1.0/6.0 * (((x1 - x2)^2)+ ((y1 - y2)^2)-2.0* ((z1 - z2)^2))*R

                    term2 = -1.0*(x1-x2)*(y1-y2)*abs(z1-z2)*atan((x1-x2)*(y1-y2),(abs(z1-z2)*R) )

                    if abs((y1-y2)+R)<1e-16
                        term3 = 0.0
                    else
                        term3 = 0.5*(((x1 - x2)^2) -((z1 - z2)^2))*(y1-y2)*real(log(complex((y1-y2)+R)))
                    end
                    if abs((x1-x2)+R)<1e-16
                        term4 = 0.0
                    else
                        term4 = 0.5*(((y1 - y2)^2) - ((z1 - z2)^2))*(x1-x2)*real(log(complex((x1-x2)+R)))
                    end

                    sol = sol + ((-1.0)^(c1+c2+c3+c4 - 4))*(term1+term2+term3+term4)
                end
            end
        end
    end

    return sol
end


function integ_surf_surf_ortho(x1v,y1v,z1,x2v,y2,z2v)
    sol=0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y1 = y1v[c2]
            for c3 in range(1, stop=2)
                z2 = z2v[c3]
                for c4 in range(1, stop=2)
                    x2 = x2v[c4]

                    R=sqrt( ((x1-x2)^2)+((y1-y2)^2)+((z1-z2)^2 ) )

                    term1 = -1.0/3.0* (y1-y2)*(z1-z2)*R

                    term2 = -1.0/2.0 * (x1-x2)*abs(z1-z2)*(z1-z2)*atan( (x1-x2)*(y1-y2),(abs(z1-z2)*R))

                    term3 = -0.5 * (x1-x2)*abs(y1-y2)*(y1-y2)*atan( (x1-x2)*(z1-z2),( abs(y1-y2)*R))

                    term4 = -1.0/6.0 * ((abs(x1-x2))^3)*atan( (y1-y2)*(z1-z2),( abs(x1-x2)*R))

                    if abs((x1-x2)+R)<1e-16
                        term5 = 0.0
                    else
                        term5 = (x1-x2)*(y1-y2)*(z1-z2)*real(log(complex( (x1-x2)+R)))
                    end

                    if abs((y1-y2)+R)<1e-16
                        term6=0.0
                    else
                        term6 = (0.5*((x1-x2)^2)-1.0/6.0*((z1-z2)^2))*(z1-z2)*real(log(complex( (y1-y2)+R)))
                    end

                    if abs((z1-z2)+R)<1e-16
                        term7=0.0
                    else
                        term7 = (0.5*((x1-x2)^2)-1.0/6.0*((y1-y2)^2))*(y1-y2)*real(log(complex( (z1-z2)+R)))
                    end

                    sol = sol + ((-1.0)^(c1+c2+c3+c4 - 4))*(term1+term2+term3+term4+term5+term6+term7)
                end
            end
        end
    end

    return sol
end


function integ_vol_surf( x1v,y1v,z1v, x2v,y2v,z2)
    sol = 0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y1 = y1v[c2]
            for c3 in range(1, stop=2)
                z1 = z1v[c3]
                for c4 in range(1, stop=2)
                    x2 = x2v[c4]
                    for c5 in range(1, stop=2)
                        y2 = y2v[c5]
                        R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                        term1 = (((z1 - z2)^2) / 12.0 - (((x1 - x2)^2) + ((y1 - y2)^ 2)) / 8.0) * (z1 - z2) * R

                        if abs((x1 - x2) + R) < 1e-16
                            term2 = 0.0
                        else
                            term2 = (((y1 - y2)^2) / 2.0 - ((z1 - z2)^2) / 6.0) * (x1 - x2) * (z1 - z2) * real(log(complex((x1-x2)+R)))
                        end

                        if abs((y1 - y2) + R) < 1e-16
                            term3 = 0.0
                        else
                            term3 = (((x1 - x2)^2) / 2.0 - ((z1 - z2)^2) / 6.0) * (y1 - y2) * (z1 - z2) * 
                                real(log(complex((y1 - y2) + R)))
                        end

                        if abs((z1 - z2) + R) < 1e-16
                            term4 = 0.0
                        else
                            term4 = (- ((x1 - x2)^4) / 24.0 - ((y1 - y2)^4) / 24.0 + 
                                 ( ((y1 - y2)^2) * ((x1 - x2)^2) ) / 4.0) *
                                real(log(complex((z1 - z2) + R)))
                        end

                        if (isnan(term4) == true || isinf(term4) == true)
                            term4 = 0.0
                        end

                        term5 = - abs(((x1 - x2)^3)) * (y1 - y2) / 6.0 * atan((y1 - y2) * (z1 - z2), (abs(x1 - x2) * R))

                        term6 = - (x1 - x2) * abs(((y1 - y2)^3)) / 6.0 * atan((x1 - x2) * (z1 - z2), (abs(y1 - y2) * R))

                        term7 = - (x1 - x2) * (y1 - y2) * abs(z1 - z2) * (z1 - z2) / 2.0 * atan((x1 - x2) * (y1 - y2),
                                                                                               (abs(z1 - z2) * R))

                        sol = sol + ((-1)^(c1 + c2 + c3 + c4 + c5 + 1 - 5)) * (
                                    term1 + term2 + term3 + term4 + term5 + term6 + term7)
                    end
                end
            end
        end
    end

    return sol
end

function integ_vol_vol(x1v,y1v,z1v, x2v,y2v,z2v)
    sol = 0.0
    for c1 in range(1, stop=2)
        x1 = x1v[c1]
        for c2 in range(1, stop=2)
            y1 = y1v[c2]
            for c3 in range(1, stop=2)
                z1 = z1v[c3]
                for c4 in range(1, stop=2)
                    x2 = x2v[c4]
                    for c5 in range(1, stop=2)
                        y2 = y2v[c5]
                        for c6 in range(1, stop=2)
                            z2 = z2v[c6]

                            R = sqrt(((x1 - x2)^2) + ((y1 - y2)^2) + ((z1 - z2)^2))

                            term1 = (((x1 - x2)^4) + ((y1 - y2)^4) + ((z1 - z2)^4)
                                     - 3.0 * ((x1 - x2)^2) * ((y1 - y2)^2) 
                                     - 3.0 * ((y1 - y2)^2) * ((z1 - z2)^2) 
                                     - 3.0 * ((x1 - x2)^2) * ((z1 - z2)^2)) * R / 60.0

                            if abs((x1 - x2) + R) < 1e-16
                                term2 = 0.0
                            else
                                term2 = (((y1 - y2)^2) * ((z1 - z2)^2) / 4.0 - ((y1 - y2)^4) / 24.0
                                     - ((z1 - z2)^4) / 24.0) *
                                    (x1 - x2) * real(log(complex((x1-x2)+R)))
                            end

                            if abs((z1 - z2) + R) < 1e-16
                                term3 = 0.0
                            else
                                term3 = (((y1 - y2)^2) * ((x1 - x2)^2) / 4.0 - ((y1 - y2)^4) / 24.0
                                    - ((x1 - x2)^4) / 24.0) * (z1 - z2) * real(log(complex((z1-z2)+R)))
                            end

                            if abs((y1 - y2) + R) < 1e-16
                                term4 = 0.0
                            else
                                term4 = (((z1 - z2)^2) * ((x1 - x2)^2) / 4.0 - ((z1 - z2)^4)  / 24.0
                                     - ((x1 - x2)^4) / 24.0) * ( y1 - y2) * real(log(complex((y1-y2)+R)))
                            end

                            term5 = - abs(((x1 - x2)^3))  * (y1 - y2) * (z1 - z2) / 6.0 * atan((y1 - y2) * (z1 - z2),(abs(x1 - x2) * R))

                            term6 = - (x1 - x2) * abs(((y1 - y2)^3)) * (z1 - z2) / 6.0 * atan((x1 - x2) * (z1 - z2),(abs(y1 - y2) * R))

                            term7 = - (x1 - x2) * (y1 - y2) * abs(((z1 - z2)^3)) / 6.0 * atan((x1 - x2) * (y1 - y2), (abs(z1 - z2) * R))
                            

                            sol = sol + ((-1.0)^(c1 + c2 + c3 + c4 + c5 + c6 + 1 - 6)) * ( term1 + term2 + term3 + term4 + term5 + term6 + term7)
                        end
                    end
                end
            end
        end
    end

    return sol
end


function compute_Lp_matrix_1(bars,sizey,sizez)
    N = size(bars)[1]
    mat_Lp=zeros((N, N))
    S = sizey * sizez * sizey * sizez
    @floop for m in range(1, stop=N)
        m1 = bars[m,1]
        m2 = bars[m,2]
        m3 = bars[m,3]
        m4 = bars[m,4]
        m5 = bars[m,5]
        m6 = bars[m,6]
        for n in range(m, stop=N)
            if m==n
                l = abs(m4 - m1)
                w = abs(m5 - m2)
                t = abs(m6 - m3)
                @inbounds mat_Lp[m, n] = compute_Lp_self(l, w, t)
            else
                @inbounds mat_Lp[m, n] = 1e-7 / S * Integ_vol_vol([m1,m4], [m2,m5], [m3,m6],
                                                    [bars[n,1],bars[n,4]], [bars[n,2],bars[n,5]], [bars[n,3],bars[n,6]])
                @inbounds mat_Lp[n, m] = mat_Lp[m, n]
            end
        end
    end
    return mat_Lp
end


function compute_Lp_matrix_2(bars,sizex,sizez)
    N = size(bars)[1]
    mat_Lp=zeros((N, N))
    S = sizex * sizez * sizex * sizez
    @floop for m in range(1, stop=N)
        m1 = bars[m,1]
        m2 = bars[m,2]
        m3 = bars[m,3]
        m4 = bars[m,4]
        m5 = bars[m,5]
        m6 = bars[m,6]
        for n in range(m, stop=N)
            if m==n
                w = abs(m4 - m1)
                l = abs(m5 - m2)
                t = abs(m6 - m3)
                @inbounds mat_Lp[m, n] = compute_Lp_self(l, w, t)
            else
                @inbounds mat_Lp[m, n] = 1e-7 / S * Integ_vol_vol([m1,m4], [m2,m5], [m3,m6],
                                                    [bars[n,1],bars[n,4]], [bars[n,2],bars[n,5]], [bars[n,3],bars[n,6]])
                @inbounds mat_Lp[n, m] = mat_Lp[m, n]
            end
        end
    end
    return mat_Lp
end


function compute_Lp_matrix_3(bars,sizex,sizey)
    N = size(bars)[1]
    mat_Lp=zeros((N, N))
    S = sizex * sizey * sizex * sizey
    @floop for m in range(1, stop=N)
        m1 = bars[m,1]
        m2 = bars[m,2]
        m3 = bars[m,3]
        m4 = bars[m,4]
        m5 = bars[m,5]
        m6 = bars[m,6]
        for n in range(m, stop=N)
            if m==n
                t = abs(m4 - m1)
                w = abs(m5 - m2)
                l = abs(m6 - m3)
                @inbounds mat_Lp[m, n] = compute_Lp_self(l, w, t)
            else
                @inbounds mat_Lp[m, n] = 1e-7 / S * Integ_vol_vol([m1,m4], [m2,m5], [m3,m6],
                                                    [bars[n,1],bars[n,4]], [bars[n,2],bars[n,5]], [bars[n,3],bars[n,6]])
                @inbounds mat_Lp[n, m] = mat_Lp[m, n]
            end
        end
    end
    return mat_Lp
end


function compute_Lps(barsLpx, barsLpy, barsLpz, sizex, sizey, sizez)
    Lp1 = compute_Lp_matrix_1(barsLpx, sizey, sizez)
    Lp2 = compute_Lp_matrix_2(barsLpy, sizex, sizez)
    Lp3 = compute_Lp_matrix_3(barsLpz, sizex, sizey)
    return Lp1, Lp2, Lp3
end

# function Lp(barsLpx, barsLpy, barsLpz, sizex, sizey, sizez, i)
#     if(i == 1)
#         return compute_Lp_matrix_1(barsLpx, sizey, sizez)
#     elseif(i == 2)
#         return compute_Lp_matrix_2(barsLpy, sizex, sizez)
#     else
#         return compute_Lp_matrix_3(barsLpz, sizex, sizey)
#     end
# end

# function compute_Lps_parallel(barsLpx, barsLpy, barsLpz, sizex, sizey, sizez)
#     addprocs(3)
#     Lps = @sync @distributed (append!) for i in 1:3
#         Lp(barsLpx, barsLpy, barsLpz, sizex, sizey, sizez, i)
#     end
#     return Lps[1], Lps[2], Lps[3]
# end

function check_condition_P(eps1,eps2,eps3,eps4,Sup1,Sup2,max_d,min_R,size_dim,other_dim1,other_dim2)
    max_oth=other_dim1*other_dim2
    condS = max(max_oth, size_dim) / (min_R + 1e-15)
    condX1a = Sup1*Sup2*max_d/((min_R+1e-15)^3)
    condX1f = size_dim/(min_R+1e-15)
    condX1b = size_dim/max_oth
    supp_dim=false
    if ((condS<eps4) && ((condX1b <= eps3 || condX1f<eps1) && condX1a<eps2))
        supp_dim=true
    end
    return supp_dim
end

function Integ_sup_sup(xc1, yc1, zc1, xc2, yc2, zc2, a1,b1,c1, a2,b2,c2)
    epsilon1 = 5e-3
    epsilon2 = 1e-3
    epsilon3 = 1e-3
    epsilon4 = 3e-1

    x1v = [xc1 - a1 / 2.0, xc1 + a1 / 2.0]
    y1v = [yc1 - b1 / 2.0, yc1 + b1 / 2.0]
    z1v = [zc1 - c1 / 2.0, zc1 + c1 / 2.0]

    x2v = [xc2 - a2 / 2.0, xc2 + a2 / 2.0]
    y2v = [yc2 - b2 / 2.0, yc2 + b2 / 2.0]
    z2v = [zc2 - c2 / 2.0, zc2 + c2 / 2.0]

    sup1_xz_plane = false
    sup1_yz_plane = false

    sup2_xz_plane = false
    sup2_yz_plane = false

    if (a1 <= b1 && a1 <= c1)
        sup1_yz_plane = true
        a1 = 1.0
    else
        if (b1 <= a1 && b1 <= c1)
            sup1_xz_plane = true
            b1 = 1.0
        else
            c1 = 1.0
        end
    end

    if (a2 <= b2 && a2 <= c2)
        sup2_yz_plane = true
        a2 = 1.0
    else
        if (b2 <= a2 && b2 <= c2)
            sup2_xz_plane = true
            b2 = 1.0
        else
            c2 = 1.0
        end
    end

    sup1 = a1 * b1 * c1
    sup2 = a2 * b2 * c2

    supp_x1 = false
    supp_y1 = false
    supp_z1 = false
    supp_x2 = false
    supp_y2 = false
    supp_z2 = false

    aux_x = [abs(x1v[1] - x2v[1]), abs(x1v[1] - x2v[2]), abs(x1v[2] - x2v[1]), abs(x1v[2] - x2v[2])]
    aux_y = [abs(y1v[1] - y2v[1]), abs(y1v[1] - y2v[2]), abs(y1v[2] - y2v[1]), abs(y1v[2] - y2v[2])]
    aux_z = [abs(z1v[1] - z2v[1]), abs(z1v[1] - z2v[2]), abs(z1v[2] - z2v[1]), abs(z1v[2] - z2v[2])]

    min_R = sqrt((minimum(aux_x)^2) + (minimum(aux_y)^2) + (minimum(aux_z)^2))

    max_d = maximum(aux_x)
    supp_x1 = check_condition_P(epsilon1, epsilon2, epsilon3, epsilon4, sup1, sup2, max_d, min_R, a1, b1, c1)
    supp_x2 = check_condition_P(epsilon1, epsilon2, epsilon3, epsilon4, sup1, sup2, max_d, min_R, a2, b2, c2)

    max_d = maximum(aux_y)
    supp_y1 = check_condition_P(epsilon1, epsilon2, epsilon3, epsilon4, sup1, sup2, max_d, min_R, b1, a1, c1)
    supp_y2 = check_condition_P(epsilon1, epsilon2, epsilon3, epsilon4, sup1, sup2, max_d, min_R, b2, a2, c2)

    max_d = maximum(aux_z)
    supp_z1 = check_condition_P(epsilon1, epsilon2, epsilon3, epsilon4, sup1, sup2, max_d, min_R, c1, a1, b1)
    supp_z2 = check_condition_P(epsilon1, epsilon2, epsilon3, epsilon4, sup1, sup2, max_d, min_R, c2, a2, b2)

    if (sup1_yz_plane == true)
        supp_x1 = true
    else
        if (sup1_xz_plane == true)
            supp_y1 = true
        else
            supp_z1 = true
        end
    end

    if (sup2_yz_plane == true)
        supp_x2 = true
    else
        if (sup2_xz_plane == true)
            supp_y2 = true
        else
            supp_z2 = true
        end
    end

    sum_supp = supp_x1 + supp_y1 + supp_z1 + supp_x2 + supp_y2 + supp_z2

    # point - point to use
    if sum_supp == 6
        integ = sup1 * sup2 / sqrt( ((xc1 - xc2)^2) + ((yc1 - yc2)^2)  + ((zc1 - zc2)^2) )
    elseif sum_supp == 5 # point - line to use
        is_point_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 3)
            is_point_v1 = true
        end

        if is_point_v1 == true
            if supp_x2 == 0 # line of volume 2 along x
                integ = sup1 * sup2 / a2 * integ_line_point(x2v, yc2, zc2, xc1, yc1, zc1)
            else
                if supp_y2 == 0 # line  of volume 2 along y
                    integ = sup1 * sup2 / b2 * integ_line_point(y2v, xc2, zc2, yc1, xc1, zc1)
                else # line of volume 2 along z
                    integ = sup1 * sup2 / c2 * integ_line_point(z2v, xc2, yc2, zc1, xc1, yc1)
                end
            end
        else
            if supp_x1 == 0 # line of volume 1 along x
                integ = sup1 * sup2 / a1 * integ_line_point(x1v, yc1, zc1, xc2, yc2, zc2)
            else
                if supp_y1 == 0 # line of volume 1 along y
                    integ = sup1 * sup2 / b1 * integ_line_point(y1v, xc1, zc1, yc2, xc2, zc2)
                else # line of volume 1 along z
                    integ = sup1 * sup2 / c1 * integ_line_point(z1v, xc1, yc1, zc2, xc2, yc2)
                end
            end
        end
    elseif sum_supp == 4  # point-surface or line-line case
        is_point_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 3)
            is_point_v1 = true
        end

        is_point_v2 = false
        if (supp_x2 + supp_y2 + supp_z2 == 3)
            is_point_v2 = true
        end
        if is_point_v1 == true #point-surface case
            if supp_x2 == true #surface of volume 2 in yz plane
                integ = sup1 * a2 * integ_point_sup(zc1, yc1, xc1, z2v, y2v, xc2)
            else
                if supp_y2 == true #surface of volume 2 in xz plane
                    integ = sup1 * b2 * integ_point_sup(xc1, zc1, yc1, x2v, z2v, yc2)
                else
                    integ = sup1 * c2 * integ_point_sup(xc1, yc1, zc1, x2v, y2v, zc2)
                end
            end
        else
            if is_point_v2 == true #point-surface case
                if supp_x1 == true #surface of volume 1 in yz plane
                    integ = sup2 * a1 * integ_point_sup(zc2, yc2, xc2, z1v, y1v, xc1)
                else
                    if supp_y1== true #surface of volume 1 in xz plane
                        integ = sup2 * b1 * integ_point_sup(xc2, zc2, yc2, x1v, z1v, yc1)
                    else #surface of volume 1 in xy plane
                        integ=sup2*c1*integ_point_sup(xc2,yc2,zc2,x1v,y1v,zc1)
                    end
                end
            else #line-line case
                if supp_y1 == true && supp_z1 == true
                    if supp_y2 == true && supp_z2 == true #parallel lines
                        integ = b1 * c1 * b2 * c2 * integ_line_line_parall(x1v, yc1, zc1,x2v, yc2, zc2)
                    else
                        if supp_x2 == true && supp_z2 == true #orthogonal lines
                            integ = b1 * c1 * a2 * c2 * integ_line_line_ortho_xy(x1v, yc1, zc1, xc2, y2v, zc2)
                        else
                            integ = b1 * c1 * a2 * b2 * integ_line_line_ortho_xy(x1v, zc1, yc1, xc2, z2v, yc2)
                        end
                    end
                else
                    if supp_x1 == true && supp_z1 == true
                        if supp_x2 == true && supp_z2 == true  # parallel lines
                            integ = a1 * c1 * a2 * c2 * integ_line_line_parall(y1v, xc1, zc1,y2v, xc2, zc2)
                        else
                            if supp_x2 == true && supp_y2 == true  # orthogonal lines
                                integ = a1 * c1 * a2 * b2 * integ_line_line_ortho_xy(y1v, zc1, xc1, yc2, z2v, xc2)
                            else
                                integ = a1 * c1 * b2 * c2 * integ_line_line_ortho_xy(y1v, xc1, zc1, yc2, x2v, zc2)
                            end
                        end
                    else
                        if supp_x2 == true && supp_y2 == true  # parallel lines
                            integ = a1 * b1 * a2 * b2 * integ_line_line_parall(z1v, xc1, yc1, z2v, xc2, yc2)
                        else
                            if supp_x2 == true && supp_z2 == true  # orthogonal lines
                                integ = a1 * b1 * a2 * c2 * integ_line_line_ortho_xy(z1v, yc1, xc1, zc2, y2v, xc2)
                            else
                                integ = a1 * b1 * b2 * c2 * integ_line_line_ortho_xy(z1v, xc1, yc1, zc2, x2v, yc2)
                            end
                        end
                    end
                end
            end
        end
    elseif sum_supp == 3 #surface-line
        is_surf_v1 = false
        if (supp_x1 + supp_y1 + supp_z1 == 1)
            is_surf_v1 = true
        end
        #line-surface case
        if is_surf_v1 == true #bar1 is a surface
            if supp_x1 == true #bar1 is a surface in y-z plane
                if supp_x2 == false #bar2 is a line along x
                    integ = a1 * b2 * c2 * integ_line_surf_ortho(x2v, yc2, zc2, xc1, y1v, z1v)
                else
                    if supp_y2 == false #bar2 is a line along y
                        integ = a1 * a2 * c2 * integ_line_surf_para(y1v, z1v, xc1, y2v, zc2, xc2)
                    else #bar2 is a line along z
                        integ = a1 * a2 * b2 * integ_line_surf_para(z1v, y1v, xc1, z2v, yc2, xc2)
                    end
                end
            else
                if supp_y1 == true #bar1 is a surface in x-z plane
                    if supp_x2 == false #bar2 is a line along x
                        integ = b1 * b2 * c2 * integ_line_surf_para(x1v, z1v, yc1, x2v, zc2, yc2)
                    else
                        if supp_y2 == false #bar2 is a line along y
                            integ = b1 * a2 * c2 * integ_line_surf_ortho(y2v, xc2, zc2, yc1, x1v, z1v)
                        else #bar2 is a line along z
                            integ = b1 * a2 * b2 * integ_line_surf_para(z1v, x1v, yc1, z2v, xc2, yc2)
                        end
                    end
                else
                    if supp_x2 == false #bar2 is a line along x
                        integ = c1 * b2 * c2 * integ_line_surf_para(x1v, y1v, zc1, x2v, yc2, zc2)
                    else
                        if supp_y2 == false #bar2 is a line along y
                            integ = c1 * a2 * c2 * integ_line_surf_para(y1v, x1v, zc1, y2v, xc2, zc2)
                        else #bar2 is a line along z
                            integ = c1 * a2 * b2 * integ_line_surf_ortho(z2v, xc2, yc2, zc1, x1v, y1v)
                        end
                    end
                end
            end
        else #bar2 is a surface
            if supp_x2 == true #bar2 is a surface in y-z plane
                if supp_x1 == false #bar1 is a line along x
                    integ=a2*b1*c1 * integ_line_surf_ortho(x1v,yc1,zc1, xc2,y2v,z2v)
                else
                    if supp_y1 == false #bar1 is a line along y
                        integ=a2*a1*c1 *integ_line_surf_para( y2v,z2v,xc2,y1v,zc1,xc1)
                    else #bar1 is a line along z
                        integ=a2*a1*b1 *integ_line_surf_para( z2v,y2v,xc2,z1v,yc1,xc1)
                    end
                end
            else
                if supp_y2 == true #bar2 is a surface in x-z plane
                    if supp_x1 == false #bar1 is a line along x
                        integ=b2*b1*c1 * integ_line_surf_para( x2v,z2v,yc2,x1v,zc1,yc1)
                    else
                        if supp_y1 == false #bar1 is a line along y
                            integ=b2*a1*c1 *integ_line_surf_ortho(y1v,xc1,zc1, yc2,x2v,z2v)
                        else #bar1 is a line along z
                            integ=b2*a1*b1 *integ_line_surf_para( z2v,x2v,yc2,z1v,xc1,yc1)
                        end
                    end
                else
                    if supp_x1 == false #bar1 is a line along x
                        integ=c2*b1*c1 * integ_line_surf_para( x2v,y2v,zc2,x1v,yc1,zc1)
                    else
                        if supp_y1 == false #bar1 is a line along y
                            integ=c2*a1*c1 *integ_line_surf_para( y2v,x2v,zc2,y1v,xc1,zc1)
                        else #bar1 is a line along z
                            integ=c2*a1*b1 *integ_line_surf_ortho(z1v,xc1,yc1, zc2,x2v,y2v)
                        end
                    end
                end
            end
        end
    else #surface-surface case
        if supp_x1 == true #bar1 is a surface in yz plane
            if supp_x2 == true #bar2 is a surface in yz plane
                integ = a1 * a2 * integ_surf_surf_para(y1v, z1v, xc1, y2v, z2v, xc2)
            else
                if supp_y2 == true #bar2 is a surface in xz plane
                    integ = a1 * b2 * integ_surf_surf_ortho(z1v, y1v, xc1, z2v, yc2, x2v)
                else #bar2 is a surface in xy plane
                    integ = a1 * c2 * integ_surf_surf_ortho(y1v, z1v, xc1, y2v, zc2, x2v)
                end
            end
        else
            if supp_y1 == true #%bar1 is a surface in xz plane
                if supp_x2 == true #bar2 is a surface in yz plane
                    integ = b1 * a2 * integ_surf_surf_ortho(z1v, x1v, yc1, z2v, xc2, y2v)
                else
                    if supp_y2 == true #bar2 is a surface in xz plane
                        integ = b1 * b2 * integ_surf_surf_para(x1v, z1v, yc1, x2v, z2v, yc2)
                    else #bar2 is a surface in xy plane
                        integ = b1 * c2 * integ_surf_surf_ortho(x1v, z1v, yc1, x2v, zc2, y2v)
                    end
                end
            else #bar1 is a surface in xy plane
                if supp_x2 == true #bar2 is a surface in yz plane
                    integ = c1 * a2 * integ_surf_surf_ortho(y1v, x1v, zc1, y2v, xc2, z2v)
                else
                    if supp_y2 == true #bar2 is a surface in xz plane
                        integ = c1 * b2 * integ_surf_surf_ortho(x1v, y1v, zc1, x2v, yc2, z2v)
                    else #bar2 is a surface in xy plane
                        integ = c1 * c2 * integ_surf_surf_para(x1v, y1v, zc1, x2v, y2v, zc2)
                    end
                end
            end
        end
    end
    res=real(integ)
    return res
end

function compute_P_matrix(centers,sup_type,sx,sy,sz)
    eps0=8.854187816997944e-12
    N=size(centers)[1]
    P_mat=zeros((N, N))
    @floop for m in range(1, stop=N)
        if sup_type[m] == 1
            A1 = sx * sz
            a1 = sx
            b1 = 0.0
            c1 = sz
        elseif sup_type[m] == 2
            A1 = sy * sz
            a1 = 0.0
            b1 = sy
            c1 = sz
        elseif sup_type[m] == 3
            A1 = sx * sy
            a1 = sx
            b1 = sy
            c1 = 0.0
        end
        for n in range(m, stop=N)
            if sup_type[n] == 1
                A2 = sx * sz
                a2 = sx
                b2 = 0.0
                c2 = sz
            elseif sup_type[n] == 2
                A2 = sy * sz
                a2 = 0.0
                b2 = sy
                c2 = sz
            elseif sup_type[n] == 3
                A2 = sx * sy
                a2 = sx
                b2 = sy
                c2 = 0.0
            end
            P_mat[m,n]=1.0 / (4.0 * pi * eps0 * A1 * A2) * Integ_sup_sup(
                centers[m,1], centers[m,2], centers[m,3], centers[n,1], centers[n,2], centers[n,3],
                a1, b1, c1, a2, b2, c2)
            P_mat[n, m]=P_mat[m,n]
        end
    end
    return P_mat
end