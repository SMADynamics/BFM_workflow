import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig, inv



def stretch_ellipse_selected_pts(x,y, p0=None, p1=None, stretchit=True, plots=False, path_save=None):
    ''' fit ellipse on selected points x[p0:p1], y[p0:p1], 
    useful when x,y have pauses 
    decide if stretch it to be circular, center to origin. 
    When stretchit=True stretches the short axis to become equal to the long axis '''
    # fit ellipse on partial xy [p0:p1]:
    xx0, yy0, center0, a, b, phi = makeBestEllipse(x[p0:p1], y[p0:p1], nel=50)
    # rotate full xy [:] so major axis is vertical: (check if this is always the case :NO!)
    x_rot, y_rot = rotateArray(np.array((x,y)), -phi)
    # fit again an ellipse on the rotated partial xy [p0:p1]:
    xx, yy, center, a, b, phi = makeBestEllipse(x_rot[p0:p1], y_rot[p0:p1], nel=50) 
    # translate to 0,0 and scale x to make a circle from the ellipse:
    if a>b:
        x_rot = (x_rot - center[0])
        if stretchit:
            y_rot = (y_rot - center[1])*a/b
        else:
            y_rot = (y_rot - center[1])
    else:
        if stretchit:
            x_rot = (x_rot - center[0])*b/a
        else:
            x_rot = (x_rot - center[0])
        y_rot = (y_rot - center[1])
    # fit again on the scaled circular x_rot y_rot data:
    xx, yy, center, a, b, phi = makeBestEllipse(x_rot[p0:p1], y_rot[p0:p1], nel=50) 
    if plots:
        plt.figure('stretch_ellipse_selected_pts', clear=True)
        plt.subplot(221)
        plt.plot(x, y, ',', ms=1, alpha=0.1)
        plt.plot(x[p0:p1], y[p0:p1], ',', ms=1, alpha=0.1)
        plt.plot(xx0, yy0,'-')
        plt.title('Original', fontsize=9)
        plt.axis('equal')
        plt.subplot(222)
        plt.plot(x_rot, y_rot, ',', ms=1, alpha=0.1)
        plt.plot(x_rot[p0:p1], y_rot[p0:p1], ',', ms=1, alpha=0.1)
        plt.plot(xx,yy,'-')
        plt.title('Rotated', fontsize=9)
        plt.axis('equal')
        plt.subplot(212)
        plt.plot(x, alpha=0.5)
        plt.plot(y, alpha=0.5)
        plt.axvspan(p0 if p0 else 0, p1 if p1 else -1, color='gray', alpha=0.2)

    return x_rot, y_rot, a, b



def stretch_ellipse(x,y, plots=False):
    ''' fit ellipse, stretch it to be circular, center to origin. 
    stretches the short axis to become equal to the long axis '''
    # fit ellipse on xy :
    xx0,yy0,center0,a,b,phi = makeBestEllipse(x,y, nel=50)
    # rotate xy so major axis is vertical: (check if this is always the case!)
    x_rot, y_rot = rotateArray(np.array((x,y)), -phi)
    # fit again an ellipse on the rotated xy:
    xx,yy,center,a,b,phi = makeBestEllipse(x_rot, y_rot, nel=50) 
    # translate to 0,0 and scale x to make a circle from the ellipse:
    if a>b:
        x_out = (x_rot - center[0])
        y_out = (y_rot - center[1])*a/b
    else:
        x_out = (x_rot - center[0])*b/a
        y_out = (y_rot - center[1])
    if plots:
        plt.figure('stretch_ellipse', clear=True)
        plt.subplot(211)
        plt.plot(x, y, ',', ms=1, alpha=0.1)
        plt.plot(xx0, yy0,'--')
        plt.title('original + elliptical fit', fontsize=9)
        plt.axis('equal')
        plt.subplot(212)
        plt.plot(x_out, y_out, ',', ms=1, alpha=0.1)
        plt.title(f'rotated and stretched to a circle', fontsize=9)
        plt.axis('equal')
        plt.tight_layout()
    return x_out, y_out #, a, b



def makeBestEllipse(x,y, nel=100):
    ''' x,y : input position data
    nel: number of pts in best ellipse
    returns: xx,yy,center,a,b,phi parameters of best ellipse'''
    a = fitEllipse(x,y)
    center = ellipse_center(a)
    phi = ellipse_angle_of_rotation(a)
    a,b = ellipse_axis_length(a)
    epts = np.arange(0, 2*np.pi, 2*np.pi/nel)
    xx = center[0] + a*np.cos(epts)*np.cos(phi) - b*np.sin(epts)*np.sin(phi)
    yy = center[1] + a*np.cos(epts)*np.sin(phi) + b*np.sin(epts)*np.cos(phi)
    return xx,yy,center,a,b,phi



def fitEllipse__ORIGINAL(x,y):
    """
    DEPRECATED: sometimes it produces complex solution instead of real. Use fitEllipse() instead.
    
    Algorithm from Fitzgibbon et al 1996, Direct Least Squares Fitting of Ellipsees.  
    Formulated in terms of Langrangian multipliers, rewritten as a generalized eigenvalue problem. """
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T, D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
    


def fitEllipse__GTP1(x, y):
    """GTP. Direct least-squares ellipse fit via the generalized eigenproblem."""
    from scipy.linalg import eig as generalized_eig
    # rm mean from x,y:
    # x = x - np.mean(x)
    # y = y - np.mean(y)
    x = np.asarray(x, dtype=float)[:, np.newaxis]
    y = np.asarray(y, dtype=float)[:, np.newaxis]
    D = np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = D.T @ D
    C = np.zeros((6, 6), dtype=float)
    C[0, 2] = 2.0
    C[2, 0] = 2.0
    C[1, 1] = -1.0
    # Solve S a = lambda C a directly, instead of eig(inv(S) @ C).
    evals, evecs = generalized_eig(S, C, check_finite=False)
    evals = np.real_if_close(evals, tol=1000)
    evecs = np.real_if_close(evecs, tol=1000)
    valid = []
    for idx in range(evecs.shape[1]):
        vec = evecs[:, idx]
        # Ignore genuinely complex solutions.
        if np.iscomplexobj(vec) and not np.allclose(np.imag(vec), 0.0, atol=1e-10):
            continue
        vec = np.real(vec)
        A, B, Cc = vec[0], vec[1], vec[2]
        # Ellipse constraint: 4AC - B^2 > 0
        if 4.0 * A * Cc - B * B > 0:
            lam = evals[idx]
            if np.iscomplexobj(lam) and not np.isclose(np.imag(lam), 0.0, atol=1e-10):
                continue
            lam = float(np.real(lam))
            if np.isfinite(lam):
                valid.append((abs(lam), vec))
    if not valid:
        raise ValueError("fitEllipse(): no valid real ellipse solution found")
    # Pick the valid ellipse mode with the smallest absolute generalized eigenvalue.
    _, a = min(valid, key=lambda item: item[0])
    # Normalize for stable downstream formulas.
    return a / np.linalg.norm(a)
    
    

def fitEllipse(x, y):
    """
    GTP. Stable direct least-squares ellipse fit (Halir-Flusser style).
    Returns conic coefficients [A, B, C, D, E, F] for:
        A x^2 + B x y + C y^2 + D x + E y + F = 0
    
    TODO if complex output, add:
            tol_imag = 1e-10
            eps_cond = 1e-12
            candidates = []
            for j in range(eigvecs.shape[1]):
                v = eigvecs[:, j]
                if np.iscomplexobj(v):
                    if np.max(np.abs(np.imag(v))) > tol_imag:
                        continue
                    v = np.real(v)
                else:
                    v = np.asarray(v, dtype=float)
                cond = 4.0 * v[0] * v[2] - v[1] * v[1]
                if cond > eps_cond and np.all(np.isfinite(v)):
                    candidates.append(v)
            if not candidates:
                raise ValueError("fitEllipse(): no valid real ellipse eigenvector")
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    if x.size != y.size or x.size < 6:
        raise ValueError("fitEllipse(): need same-length x,y with at least 6 points")

    # Normalize for conditioning
    mx, my = np.mean(x), np.mean(y)
    sx, sy = np.std(x), np.std(y)
    if sx == 0 or sy == 0:
        raise ValueError("fitEllipse(): degenerate input (zero variance)")

    xn = (x - mx) / sx
    yn = (y - my) / sy

    D1 = np.column_stack((xn * xn, xn * yn, yn * yn))
    D2 = np.column_stack((xn, yn, np.ones_like(xn)))

    S1 = D1.T @ D1
    S2 = D1.T @ D2
    S3 = D2.T @ D2

    # T = -inv(S3) S2^T, M = inv(C1) (S1 + S2 T)
    T = -np.linalg.solve(S3, S2.T)
    M = S1 + S2 @ T
    C1 = np.array([[0.0, 0.0, 2.0],
                   [0.0, -1.0, 0.0],
                   [2.0, 0.0, 0.0]])

    eigvals, eigvecs = np.linalg.eig(np.linalg.solve(C1, M))
    eigvecs = np.real_if_close(eigvecs, tol=1000)

    # Keep only ellipse solutions: 4ac - b^2 > 0
    cond = 4.0 * eigvecs[0, :] * eigvecs[2, :] - eigvecs[1, :] ** 2
    idx = np.where(cond > 0)[0]
    if idx.size == 0:
        # Soft fallback for near-degenerate numerics
        idx = [int(np.argmax(cond))]
        if cond[idx[0]] <= 0:
            raise ValueError("fitEllipse(): no valid real ellipse solution found")

    a1 = np.real(eigvecs[:, idx[0]])
    a2 = T @ a1
    an = np.concatenate((a1, a2))  # normalized-space conic coeffs

    # Denormalize coefficients back to original x,y
    A, B, C, D, E, F = an
    Ao = A / (sx * sx)
    Bo = B / (sx * sy)
    Co = C / (sy * sy)
    Do = D / sx - 2.0 * A * mx / (sx * sx) - B * my / (sx * sy)
    Eo = E / sy - 2.0 * C * my / (sy * sy) - B * mx / (sx * sy)
    Fo = (
        F
        + A * mx * mx / (sx * sx)
        + B * mx * my / (sx * sy)
        + C * my * my / (sy * sy)
        - D * mx / sx
        - E * my / sy
    )

    a = np.array([Ao, Bo, Co, Do, Eo, Fo], dtype=float)
    return a / np.linalg.norm(a)



def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f + c*d*d + g*b*b - 2*b*d*f - a*c*g)
    down1 = (b*b-a*c)*((c-a)*np.sqrt(1 + 4*b*b/((a-c)*(a-c))) - (c+a)) # getting warning from (a-c)=0
    down2 = (b*b-a*c)*((a-c)*np.sqrt(1 + 4*b*b/((a-c)*(a-c))) - (c+a))
    res1=np.sqrt(abs(up/down1))
    res2=np.sqrt(abs(up/down2))
    return np.array([res1, res2])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c)) # getting warning from (a-c)=0


def rotateArray(a,th):
    ''' rotates the array a of the angle theta (rad)'''
    R = np.array(([np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]))
    rota = np.dot(R,a)
    return rota
