function dudx=DudxForward2(uf2,uf1,ui,dx)

    dudx=(-uf2+4*uf1-3*ui)/(2*dx);

end