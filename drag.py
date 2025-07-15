import numpy as np

eta = 0.001   # [Pa*s] = [N*s/m^2] water viscosity



def cal_drag(traj_radius_m, bead_diam_m, dist_beadsurf_wall):
    ''' calculate corrected drag of a rotating and translating bead (moon-like)
            - of radius "bead_radius" [m], 
            - rotating on a horizontal circular trajectory whose center is displaced by "axis_offset" [m] from the bead center,
            - at a vertical distance from the wall (bead surface - wall surface) of "dist_beadsurf_wall" [m]
        taking the maximum between Faxen and Brenner proximity corrections 
    '''
    bead_radius_m = bead_diam_m/2
    # drag [pN nm s], Faxen and Brenner corrections:
    gamma_faxen_pNnms = calc_drag_faxen(bead_radius=bead_radius_m, axis_offset=traj_radius_m, dist_beadsurf_wall=dist_beadsurf_wall)
    gamma_brenn_pNnms = calc_drag_brenn(bead_radius=bead_radius_m, axis_offset=traj_radius_m, dist_beadsurf_wall=dist_beadsurf_wall)
    # take drag as the max between Faxen and Brenner corrections:
    gamma_pNnms = np.max((gamma_brenn_pNnms, gamma_faxen_pNnms), axis=0)
    return gamma_pNnms



def calc_drag_faxen(bead_radius, axis_offset, dist_beadsurf_wall, prints=0):
    ''' returns the FAXEN drag (in pN nm s) 
        see: 'Comparison of Faxen s correction for a microsphere translating or rotating near a surface' PRE (2009)
    '''
    dist_beadcent_wall = bead_radius + dist_beadsurf_wall   # [m] dist(bead center, wall)
    # correction beta_parallel (torque)
    faxen_1 = 1 - (1/8.)*(bead_radius/dist_beadcent_wall)**3 
    # Faxen correction gamma_parallel (force):
    faxen_2   = 1 - (9/16.)*(bead_radius/dist_beadcent_wall) + (1./8)*(bead_radius/dist_beadcent_wall)**3
    faxen_2_1 = faxen_2 - (45/256.)*(bead_radius/dist_beadcent_wall)**4 - (1./16)*(bead_radius/dist_beadcent_wall)**5
    # uncorrected rotational and translational (valid for bulk volume):
    rot_bulk = 8*np.pi*eta*bead_radius**3
    transl_bulk = 6*np.pi*eta*bead_radius*axis_offset**2 
    # proximity corrected rotational and translational:
    rot = 8*np.pi*eta*bead_radius**3/faxen_1 
    trans = transl_bulk/faxen_2_1
    # full drag:
    drag_bead = rot + trans 
    drag_bead_pNnms = drag_bead*1e21
    if prints:	    
        k_hook = 400 # [pN*nm/rad] angular stiffness of the flagellar hook (see S.Block polyhook paper)
        print('using faxen_2 correction')
        print("bead drag : "+str(drag_bead_pNnms)+" pN nm s")
        print("charact.time on hook: "+str(1000*drag_bead_pNnms/k_hook)+" ms")
        print("charact.freq. on hook: "+str(1./(drag_bead_pNnms/k_hook))+" Hz")
        print('faxen 1 = '+str(faxen_1))
        print('faxen 2_1 = '+str(faxen_2_1))
        print('faxen 2 = '+str(faxen_2))
        print(f'bulk rot drag = {rot_bulk *1e21}')
        print(f'bulk trasl drag = {transl_bulk*1e21}')
        print(f'bulk trasl + rot  = {(transl_bulk + rot_bulk)*1e21}')
    return drag_bead_pNnms



def calc_drag_brenn(bead_radius, axis_offset, dist_beadsurf_wall):
    ''' Calculate drag [pN nm s] with BRENNER corrections
        see Brenner 1967 "Slow viscous motion of a sphere parallel to a plane wall"
    '''
    z3 = 1.20205 # Riemann Z(3)
    # uncorrected rotational and translational (valid for bulk volume):
    rot_bulk = 8*np.pi*eta*bead_radius**3
    trans_bulk = 6*np.pi*eta*axis_offset**2 * bead_radius
    # proximity corrected rotational and translational:
    rot = rot_bulk * (z3 - 3*(np.pi**2/6 - 1)*(dist_beadsurf_wall/bead_radius)) 
    trans = np.abs(trans_bulk * (8/15 * np.log(dist_beadsurf_wall/bead_radius) - 0.9588))
    # full drag:
    gamma = rot + trans
    gamma_pNnms = gamma*1e21
    return gamma_pNnms

