function EvalAnalDani( x, y, rc, ηm, ηc )
    # ---------------------------------------------------------------------------
    # ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION:
    #
    # BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
    # FAR FIELD FLOW - VISCOSITIES - GEOMETRY
    #
    # ---------------------------------------------------------------------------
    
    # Input:
    gr  = 0                        # Simple shear: gr=1, er=0
    er  = -1                       # Strain rate
    A   = ηm*(ηc-ηm)/(ηc+ηm)
    i   = 1im

    if sqrt(x^2.0 + y^2.0)<=rc
        # Inclusion
        Z       =   x + i*y;
        # Velocity
        V_tot      = (ηm/(ηc+ηm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z
        vx         = real(V_tot)
        vy         = imag(V_tot)
        # Pressure 
        η          = ηc
        P          = 0.
        # Stresses
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
        # Matrix
        Z            = x + i*y;
        # Pressure
        P            = -2.0*ηm.*(ηc-ηm)./(ηc+ηm).*real(rc^2.0/Z.^2.0*(i*gr+2*er));
        # Velocity
        phi_z        = -(i/2)*ηm*gr*Z-(i*gr+2*er)*A*rc^2*Z^(-1);
        d_phi_z      = -(i/2)*ηm*gr + (i*gr+2*er)*A*rc^2/Z^2;
        conj_d_phi_z = conj(d_phi_z);
        psi_z        = (i*gr-2*er)*ηm*Z-(i*gr+2*er)*A*rc^4*Z^(-3);
        conj_psi_z   = conj(psi_z);
        V_tot        = (phi_z- Z*conj_d_phi_z - conj_psi_z) / (2*ηm);
        vx           = real(V_tot);
        vy           = imag(V_tot);
        η            = ηm;
        # Stresses
        ∂vx∂x        = (2.0*er*ηm)/(ηm + ηc)
        ∂vx∂y        = (gr*(3.0*ηm + ηc))/(2.0*(ηm + ηc))
        ∂vy∂y        = -(2.0*er*ηm)/(ηm + ηc)
        ∂vy∂x        = (gr*(ηm - ηc))/(2.0*(ηm + ηc))
        τxy          =  η*(∂vx∂y + ∂vy∂x)
        τxx          = 2η* ∂vx∂x
        τyy          = 2η* ∂vy∂y
        σxx          = -P + τxx
        σyy          = -P + τyy
        σxy          = τxy
    end

    return vx, vy, P, σxx, σyy, σxy, ∂vx∂x, ∂vx∂y, ∂vy∂x, ∂vy∂y, η
end

#--------------------------------------------------------------------#

function Tractions( x, y, rc, ηm, ηc, phase )

    # Input:
    gr  = 0;                        # Simple shear: gr=1, er=0
    er  = -1;                       # Strain rate
    A   =   ηm*(ηc-ηm)/(ηc+ηm);
    i   =   1im;

    if phase == 2
        # Inclusion
        Z          =   x + i*y;
        # Velocity
        V_tot      =  (ηm/(ηc+ηm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z;
        vx         =  real(V_tot);
        vy         =  imag(V_tot);
        # Pressure        
        P          =   0.;   # if you want you can add NaN, but according to Schmid's thesis it's zero inside
        # Stresses
        ∂Vx∂x      = (2*er*ηm)/(ηm + ηc);
        ∂Vx∂y      = (gr*(3*ηm + ηc))/(2*(ηm + ηc));
        ∂Vy∂x      = (gr*(ηm - ηc))/(2*(ηm + ηc));
        ∂Vy∂y      = -(2*er*ηm)/(ηm + ηc);
        τxy        = ηm*(∂Vx∂y + ∂Vy∂x);
        τxx        = ηm*(∂Vx∂x + ∂Vx∂x);
        τyy        = ηm*(∂Vy∂y + ∂Vy∂y);
        sxx        = real(-P + τxx)
        syy        = real(-P + τyy)
        sxy        = real(τxy)
    end
    if phase == 1
        # Matrix
        Z            =   x + i*y;
        # Pressure
        P            = -2.0*ηm.*(ηc-ηm)./(ηc+ηm).*real(rc^2.0/Z.^2.0*(i*gr+2*er));
        # Velocity
        phi_z        = -(i/2)*ηm*gr*Z-(i*gr+2*er)*A*rc^2*Z^(-1);
        d_phi_z      = -(i/2)*ηm*gr + (i*gr+2*er)*A*rc^2/Z^2;
        conj_d_phi_z = conj(d_phi_z);
        psi_z        = (i*gr-2*er)*ηm*Z-(i*gr+2*er)*A*rc^4*Z^(-3);
        conj_psi_z   = conj(psi_z);
        V_tot        = (phi_z- Z*conj_d_phi_z - conj_psi_z) / (2*ηm);
        vx           = real(V_tot);
        vy           = imag(V_tot);
        # Evaluate stresses (Valid one)
        A            = ηm*(ηc-ηm)/(ηc+ηm);
        ∂Vx∂x        = (- 3*A*er*rc^4*x^4 - 6*A*gr*rc^4*x^3*y + 18*A*er*rc^4*x^2*y^2 + 6*A*gr*rc^4*x*y^3 - 3*A*er*rc^4*y^4 + 2*A*er*rc^2*x^6 + 4*A*gr*rc^2*x^5*y - 10*A*er*rc^2*x^4*y^2 - 10*A*er*rc^2*x^2*y^4 - 4*A*gr*rc^2*x*y^5 + 2*A*er*rc^2*y^6 + er*ηm*x^8 + 4*er*ηm*x^6*y^2 + 6*er*ηm*x^4*y^4 + 4*er*ηm*x^2*y^6 + er*ηm*y^8)/(ηm*(x^2 + y^2)^4);
        ∂Vx∂y        = (3*A*gr*rc^4*x^4 - 24*A*er*rc^4*x^3*y - 18*A*gr*rc^4*x^2*y^2 + 24*A*er*rc^4*x*y^3 + 3*A*gr*rc^4*y^4 - 4*A*gr*rc^2*x^6 + 24*A*er*rc^2*x^5*y + 8*A*gr*rc^2*x^4*y^2 + 16*A*er*rc^2*x^3*y^3 + 12*A*gr*rc^2*x^2*y^4 - 8*A*er*rc^2*x*y^5 + 2*gr*ηm*x^8 + 8*gr*ηm*x^6*y^2 + 12*gr*ηm*x^4*y^4 + 8*gr*ηm*x^2*y^6 + 2*gr*ηm*y^8)/(2*ηm*(x^2 + y^2)^4);
        ∂Vy∂x        = (A*rc^2*(3*gr*rc^2*x^4 - 24*er*rc^2*x^3*y - 18*gr*rc^2*x^2*y^2 + 24*er*rc^2*x*y^3 + 3*gr*rc^2*y^4 + 8*er*x^5*y + 12*gr*x^4*y^2 - 16*er*x^3*y^3 + 8*gr*x^2*y^4 - 24*er*x*y^5 - 4*gr*y^6))/(2*ηm*(x^2 + y^2)^4);
        ∂Vy∂y        = -(- 3*A*er*rc^4*x^4 - 6*A*gr*rc^4*x^3*y + 18*A*er*rc^4*x^2*y^2 + 6*A*gr*rc^4*x*y^3 - 3*A*er*rc^4*y^4 + 2*A*er*rc^2*x^6 + 4*A*gr*rc^2*x^5*y - 10*A*er*rc^2*x^4*y^2 - 10*A*er*rc^2*x^2*y^4 - 4*A*gr*rc^2*x*y^5 + 2*A*er*rc^2*y^6 + er*ηm*x^8 + 4*er*ηm*x^6*y^2 + 6*er*ηm*x^4*y^4 + 4*er*ηm*x^2*y^6 + er*ηm*y^8)/(ηm*(x^2 + y^2)^4);
        τxy          = ηc*(∂Vx∂y + ∂Vy∂x);
        τxx          = ηc*(∂Vx∂x + ∂Vx∂x);
        τyy          = ηc*(∂Vy∂y + ∂Vy∂y);
        sxx          = -P + τxx
        syy          = -P + τyy
        sxy          = τxy
    end

    return P, real(∂Vx∂x), real(∂Vx∂y), real(∂Vy∂x), real(∂Vy∂y)
end
