function dudx=DudxBackward2(ub2,ub1,ui,dx)

    dudx=(ub2-4*ub1+3*ui)/(2*dx);

end