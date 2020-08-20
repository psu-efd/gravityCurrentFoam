//solve equations for the salinity and sediment (separately)
{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    //For the suspended sediment
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh
        ),
        phi + (mesh.Sf() & fvc::interpolate(Vs))
    );


    fvScalarMatrix SedEqn
    (
        fvm::ddt(Sed)
      + fvm::div(phiAlpha, Sed)
      - fvm::laplacian(kappaEff, Sed)
    );

    SedEqn.relax();
    SedEqn.solve();

    //for the salinity
    fvScalarMatrix SalEqn
    (
        fvm::ddt(Sal)
      + fvm::div(phi, Sal)
      - fvm::laplacian(kappaEff, Sal)
    );

    SalEqn.relax();
    SalEqn.solve();

    bound(Sal,0.0);
    bound(Sed,0.0);

    rhok = calc_density(Sal,Sed);
}
