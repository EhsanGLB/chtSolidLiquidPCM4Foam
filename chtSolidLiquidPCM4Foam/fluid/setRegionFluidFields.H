    const fvMesh& mesh = fluidRegions[i];

    volScalarField& T = TFluid[i];
    volScalarField& p = pFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& Kappa = KappaFluid[i];
    incompressible::turbulenceModel& turb = turbulence[i];
    word& type = typeFluid[i];
    dimensionedScalar& rho = rhoFluid[i];
    dimensionedScalar& Cp = CpFluid[i];
    dimensionedScalar& Pr = PrFluid[i];
    dimensionedScalar& Prt = PrtFluid[i];

    mesh.schemesDict().setFluxRequired(p.name());
