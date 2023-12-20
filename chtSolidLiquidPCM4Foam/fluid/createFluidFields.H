    // Initialise fluid field pointer lists
    PtrList<volScalarField> TFluid(fluidRegions.size());
    PtrList<volScalarField> pFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<volScalarField> KappaFluid(fluidRegions.size());
    PtrList<singlePhaseTransportModel> laminarTransport(fluidRegions.size());
    PtrList<incompressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<word> typeFluid(fluidRegions.size());

    PtrList<dimensionedScalar> rhoFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CpFluid(fluidRegions.size());
    PtrList<dimensionedScalar> PrFluid(fluidRegions.size());
    PtrList<dimensionedScalar> PrtFluid(fluidRegions.size());


    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to TFluid\n" << endl;
        TFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "T",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to pFluid\n" << endl;
        pFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(UFluid[i]) & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to KappaFluid\n" << endl;
        KappaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Kappa",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Kappa", dimThermalConductivity, SMALL)
            )
        );

        Info<< "    Adding to laminarTransport\n" << endl;
        laminarTransport.set
        (
            i,
            new singlePhaseTransportModel
            (
                UFluid[i],
                phiFluid[i]
            )
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<incompressible::turbulenceModel>
            (
                incompressible::turbulenceModel::New
                (
                    UFluid[i],
                    phiFluid[i],
                    laminarTransport[i]
                )
            ).ptr()
        );

        typeFluid.set
        (
            i,
            new word (laminarTransport[i].lookup("type"))
        );

        rhoFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("rho"))
        );

        CpFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Cp"))
        );

        PrFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Pr"))
        );

        PrtFluid.set
        (
            i,
            new dimensionedScalar (laminarTransport[i].lookup("Prt"))
        );
    }




// conditions

    // forced convection
    PtrList<volScalarField> rAUFluid(fluidRegions.size());
    PtrList<pimpleControl> pimpleFluid(fluidRegions.size());


    // natural convection and PCM
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<dimensionedScalar> betaFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TRefFluid(fluidRegions.size());
    PtrList<volScalarField> rhokFluid(fluidRegions.size());


    // PCM
    PtrList<volScalarField> gammaFluid(fluidRegions.size());
    PtrList<dimensionedScalar> hmFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmFluid(fluidRegions.size());
    PtrList<dimensionedScalar> TmrFluid(fluidRegions.size());
    PtrList<dimensionedScalar> CuFluid(fluidRegions.size());
    PtrList<dimensionedScalar> bFluid(fluidRegions.size());


    forAll(fluidRegions, i)
    {
      // forced convection
      if(typeFluid[i] == "forcedConvection")
      {
        Info<< "Reading field rAU if present\n" << endl;
        rAUFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rAU",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                runTime.deltaT(),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        Info<< "    Adding to pimpleFluid\n" << endl;
        pimpleFluid.set
        (
            i,
            new pimpleControl
            (
                fluidRegions[i],
                "PIMPLE"
            )
        );
      }


      // natural convection and PCM
      if(typeFluid[i] == "naturalConvection" || "PCM")
      {
        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        if(laminarTransport[i].found("beta"))
        {
            betaFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("beta"))
            );
        }

        if(laminarTransport[i].found("TRef"))
        {
            TRefFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("TRef"))
            );
        }

        Info<< "    Adding to rhokFluid\n" << endl;
        rhokFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhok",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("rhok", dimless, SMALL)
            )
        );
      }



      // PCM
      if(typeFluid[i] == "PCM")
      {
        Info<< "    Adding to gammaFluid\n" << endl;
        gammaFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "gamma",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("gamma", dimless, SMALL)
            )
        );

        if(laminarTransport[i].found("hm"))
        {
            hmFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("hm"))
            );
        }

        if(laminarTransport[i].found("Tm"))
        {
            TmFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("Tm"))
            );
        }

        if(laminarTransport[i].found("Tmr"))
        {
            TmrFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("Tmr"))
            );
        }

        if(laminarTransport[i].found("Cu"))
        {
            CuFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("Cu"))
            );
        }

        if(laminarTransport[i].found("b"))
        {
            bFluid.set
            (
                i,
                new dimensionedScalar (laminarTransport[i].lookup("b"))
            );
        }
      }
    }
