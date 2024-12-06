function xdot = CL_Luenberger(t,x,r,AL,Acl,B,K,F)
    
    e = x(7:12);
    statedot = Acl*x(1:6) + B*K*e + B*F*r;
    edot = AL*e;

    xdot = [statedot;edot];

end

