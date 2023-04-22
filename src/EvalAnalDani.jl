function EvalAnalDani( x, y, rc, ηm, ηc )
# ---------------------------------------------------------------------------
# ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION:
#
# BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
# FAR FIELD FLOW - VISCOSITIES - GEOMETRY
#
# ---------------------------------------------------------------------------
# INPUT:
gr  = 0                        # Simple shear: gr=1, er=0
er  = -1                       # Strain rate
A   = ηm*(ηc-ηm)/(ηc+ηm)
i   = 1im
# --------------------------------------------------------------
# PRESSURE CALCULATION OUTSIDE OF AN INCLUSION IN THE y-PLANE
# --------------------------------------------------------------
# INSIDE CLAST
if sqrt(x^2.0 + y^2.0)<=rc
    Z       =   x + i*y
    P       =   0;  
    # VELOCITY
    V_tot      = (ηm/(ηc+ηm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z
    vx         = real(V_tot)
    vy         = imag(V_tot)
    η          = ηc
    A          = ηm*(ηc-ηm)/(ηc+ηm)
    ∂vx∂x      = (- 3.0*A*er*rc^4*x^4 - 6.0*A*gr*rc^4*x^3.0*y + 18.0*A*er*rc^4*x^2.0*y^2.0 + 6.0*A*gr*rc^4*x*y^3.0 - 3.0*A*er*rc^4*y^4 + 2.0*A*er*rc^2.0*x^6.0 + 4*A*gr*rc^2.0*x^5*y - 10*A*er*rc^2.0*x^4*y^2.0 - 10*A*er*rc^2.0*x^2.0*y^4 - 4*A*gr*rc^2.0*x*y^5 + 2.0*A*er*rc^2.0*y^6.0 + er*ηm*x^8 + 4*er*ηm*x^6.0*y^2.0 + 6.0*er*ηm*x^4*y^4 + 4*er*ηm*x^2.0*y^6.0 + er*ηm*y^8)/(ηm*(x^2.0 + y^2.0)^4)
    ∂vx∂y      = (3.0*A*gr*rc^4*x^4 - 24.0*A*er*rc^4*x^3.0*y - 18.0*A*gr*rc^4*x^2.0*y^2.0 + 24.0*A*er*rc^4*x*y^3.0 + 3.0*A*gr*rc^4*y^4 - 4*A*gr*rc^2.0*x^6.0 + 24.0*A*er*rc^2.0*x^5*y + 8*A*gr*rc^2.0*x^4*y^2.0 + 16*A*er*rc^2.0*x^3.0*y^3.0 + 12*A*gr*rc^2.0*x^2.0*y^4 - 8*A*er*rc^2.0*x*y^5 + 2.0*gr*ηm*x^8 + 8*gr*ηm*x^6.0*y^2.0 + 12*gr*ηm*x^4*y^4 + 8*gr*ηm*x^2.0*y^6.0 + 2.0*gr*ηm*y^8)/(2.0*ηm*(x^2.0 + y^2.0)^4)
    ∂vy∂x      = (A*rc^2.0*(3.0*gr*rc^2.0*x^4 - 24.0*er*rc^2.0*x^3.0*y - 18.0*gr*rc^2.0*x^2.0*y^2.0 + 24.0*er*rc^2.0*x*y^3.0 + 3.0*gr*rc^2.0*y^4 + 8*er*x^5*y + 12*gr*x^4*y^2.0 - 16*er*x^3.0*y^3.0 + 8*gr*x^2.0*y^4 - 24.0*er*x*y^5 - 4*gr*y^6.0))/(2.0*ηm*(x^2.0 + y^2.0)^4)
    ∂vy∂y      = -(- 3.0*A*er*rc^4*x^4 - 6.0*A*gr*rc^4*x^3.0*y + 18.0*A*er*rc^4*x^2.0*y^2.0 + 6.0*A*gr*rc^4*x*y^3.0 - 3.0*A*er*rc^4*y^4 + 2.0*A*er*rc^2.0*x^6.0 + 4*A*gr*rc^2.0*x^5*y - 10*A*er*rc^2.0*x^4*y^2.0 - 10*A*er*rc^2.0*x^2.0*y^4 - 4*A*gr*rc^2.0*x*y^5 + 2.0*A*er*rc^2.0*y^6.0 + er*ηm*x^8 + 4*er*ηm*x^6.0*y^2.0 + 6.0*er*ηm*x^4*y^4 + 4*er*ηm*x^2.0*y^6.0 + er*ηm*y^8)/(ηm*(x^2.0 + y^2.0)^4)
    τxy        =  η*(∂vx∂y + ∂vy∂x)
    τxx        = 2η* ∂vx∂x
    τyy        = 2η* ∂vy∂y
    σxx        = real(-P + τxx)
    σyy        = real(-P + τyy)
    σxy        = real(τxy)
else
    # OUTSIDE CLAST, RESP. MATRIX
    Z              =   x + i*y;
    # PRESSURE
    P            =   -2.0*ηm.*(ηc-ηm)./(ηc+ηm).*real(rc^2.0/Z.^2.0*(i*gr+2*er));
    # VELOCITY
    phi_z        = -(i/2)*ηm*gr*Z-(i*gr+2*er)*A*rc^2*Z^(-1);
    d_phi_z      = -(i/2)*ηm*gr + (i*gr+2*er)*A*rc^2/Z^2;
    conj_d_phi_z = conj(d_phi_z);
    psi_z        = (i*gr-2*er)*ηm*Z-(i*gr+2*er)*A*rc^4*Z^(-3);
    conj_psi_z   = conj(psi_z);
    V_tot        = (phi_z- Z*conj_d_phi_z - conj_psi_z) / (2*ηm);
    vx           =  real(V_tot);
    vy           =  imag(V_tot);
    η            = ηm;
    # Evaluate stress
    ∂vx∂x = (2.0*er*ηm)/(ηm + ηc)
    ∂vx∂y = (gr*(3.0*ηm + ηc))/(2.0*(ηm + ηc))
    ∂vy∂y = -(2.0*er*ηm)/(ηm + ηc)
    ∂vy∂x = (gr*(ηm - ηc))/(2.0*(ηm + ηc))
    τxy   =  η*(∂vx∂y + ∂vy∂x)
    τxx   = 2η* ∂vx∂x
    τyy   = 2η* ∂vy∂y
    σxx   = -P + τxx
    σyy   = -P + τyy
    σxy   = τxy
end

return vx, vy, P, σxx, σyy, σxy, ∂vx∂x, ∂vx∂y, ∂vy∂x, ∂vy∂y, η
end


function Tractions( x, y, rc, ηm, ηc, phase )

# ---------------------------------------------------------------------------
# ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION:
#
# BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
# FAR FIELD FLOW - VISCOSITIES - GEOMETRY
#
# ---------------------------------------------------------------------------

# INPUT:
gr  = 0;                        # Simple shear: gr=1, er=0
er  = -1;                       # Strain rate
A   =   ηm*(ηc-ηm)/(ηc+ηm);
i   =   1im;

# --------------------------------------------------------------
# PRESSURE CALCULATION OUTSIDE OF AN INCLUSION IN THE y-PLANE
# --------------------------------------------------------------

# INSIDE CLAST
if phase == 2
    
    Z       =   x + i*y;
    P       =   0;   # if you want you can add NaN, but according to Schmid's thesis it's zero inside
    
    # VELOCITY
    V_tot      =  (ηm/(ηc+ηm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z;
    vx         =  real(V_tot);
    vy         =  imag(V_tot);
    eta        =  ηc;
    
    # Evaluate stresses
    v1x = (2*er*ηm)/(ηm + ηc);
    v1y = (gr*(3*ηm + ηc))/(2*(ηm + ηc));
    v2x = (gr*(ηm - ηc))/(2*(ηm + ηc));
    v2y = -(2*er*ηm)/(ηm + ηc);

    Txy = ηm*(v1y + v2x);
    Txx = ηm*(v1x + v1x);
    Tyy = ηm*(v2y + v2y);
    sxx = real(-P + Txx)
    syy = real(-P + Tyy)
    sxy = real(Txy)
end
if phase == 1
    # OUTSIDE CLAST, RESP. MATRIX
    Z              =   x + i*y;
    # PRESSURE
    P          =   -2.0*ηm.*(ηc-ηm)./(ηc+ηm).*real(rc^2.0/Z.^2.0*(i*gr+2*er));
    
    # VELOCITY
    phi_z          = -(i/2)*ηm*gr*Z-(i*gr+2*er)*A*rc^2*Z^(-1);
    d_phi_z        = -(i/2)*ηm*gr + (i*gr+2*er)*A*rc^2/Z^2;
    conj_d_phi_z   = conj(d_phi_z);
    psi_z          = (i*gr-2*er)*ηm*Z-(i*gr+2*er)*A*rc^4*Z^(-3);
    conj_psi_z     = conj(psi_z);
    
    V_tot          = (phi_z- Z*conj_d_phi_z - conj_psi_z) / (2*ηm);
    vx         =  real(V_tot);
    vy         =  imag(V_tot);
    eta        = ηm;
      
    # Evaluate stresses (Valid one)
    A = ηm*(ηc-ηm)/(ηc+ηm);
    v1x = (- 3*A*er*rc^4*x^4 - 6*A*gr*rc^4*x^3*y + 18*A*er*rc^4*x^2*y^2 + 6*A*gr*rc^4*x*y^3 - 3*A*er*rc^4*y^4 + 2*A*er*rc^2*x^6 + 4*A*gr*rc^2*x^5*y - 10*A*er*rc^2*x^4*y^2 - 10*A*er*rc^2*x^2*y^4 - 4*A*gr*rc^2*x*y^5 + 2*A*er*rc^2*y^6 + er*ηm*x^8 + 4*er*ηm*x^6*y^2 + 6*er*ηm*x^4*y^4 + 4*er*ηm*x^2*y^6 + er*ηm*y^8)/(ηm*(x^2 + y^2)^4);
    v1y = (3*A*gr*rc^4*x^4 - 24*A*er*rc^4*x^3*y - 18*A*gr*rc^4*x^2*y^2 + 24*A*er*rc^4*x*y^3 + 3*A*gr*rc^4*y^4 - 4*A*gr*rc^2*x^6 + 24*A*er*rc^2*x^5*y + 8*A*gr*rc^2*x^4*y^2 + 16*A*er*rc^2*x^3*y^3 + 12*A*gr*rc^2*x^2*y^4 - 8*A*er*rc^2*x*y^5 + 2*gr*ηm*x^8 + 8*gr*ηm*x^6*y^2 + 12*gr*ηm*x^4*y^4 + 8*gr*ηm*x^2*y^6 + 2*gr*ηm*y^8)/(2*ηm*(x^2 + y^2)^4);
    v2x = (A*rc^2*(3*gr*rc^2*x^4 - 24*er*rc^2*x^3*y - 18*gr*rc^2*x^2*y^2 + 24*er*rc^2*x*y^3 + 3*gr*rc^2*y^4 + 8*er*x^5*y + 12*gr*x^4*y^2 - 16*er*x^3*y^3 + 8*gr*x^2*y^4 - 24*er*x*y^5 - 4*gr*y^6))/(2*ηm*(x^2 + y^2)^4);
    v2y = -(- 3*A*er*rc^4*x^4 - 6*A*gr*rc^4*x^3*y + 18*A*er*rc^4*x^2*y^2 + 6*A*gr*rc^4*x*y^3 - 3*A*er*rc^4*y^4 + 2*A*er*rc^2*x^6 + 4*A*gr*rc^2*x^5*y - 10*A*er*rc^2*x^4*y^2 - 10*A*er*rc^2*x^2*y^4 - 4*A*gr*rc^2*x*y^5 + 2*A*er*rc^2*y^6 + er*ηm*x^8 + 4*er*ηm*x^6*y^2 + 6*er*ηm*x^4*y^4 + 4*er*ηm*x^2*y^6 + er*ηm*y^8)/(ηm*(x^2 + y^2)^4);
    Txy = ηc*(v1y + v2x);
    Txx = ηc*(v1x + v1x);
    Tyy = ηc*(v2y + v2y);
    sxx = -P + Txx
    syy = -P + Tyy
    sxy = Txy
end

return P, real(v1x), real(v1y), real(v2x), real(v2y)
end
