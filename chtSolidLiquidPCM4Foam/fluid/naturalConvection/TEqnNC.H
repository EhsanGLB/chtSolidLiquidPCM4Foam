//for (int i=0; i<10; i++ ) // added
{
    volScalarField kappaEff( "kappaEff", turb.nu()/Pr + turb.nut()/Prt );
    Kappa = kappaEff*rho*Cp;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();
    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}
