for (int j=0; j<7; j++ ) // added
{
    volScalarField gammaT ("gammaT", hm*4.0*exp(-pow(4.0*(T-Tm)/Tmr,2))/Foam::sqrt(mathematicalConstant::pi)/Tmr); //added
    volScalarField kappaEff( "kappaEff", turb.nu()/Pr + turb.nut()/Prt );
    Kappa = kappaEff*rho*Cp;

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
      + (1/Cp)*gammaT*fvm::ddt(T)
      + (1/Cp)*gammaT*(U & fvc::grad(T))
    );

    TEqn.relax();
    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
    gamma = 0.5 * erf( 4.0 * ( T - Tm ) / Tmr ) + scalar(0.5);
}
