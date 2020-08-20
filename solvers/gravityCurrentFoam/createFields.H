    Info<< "Reading salinity field Sal\n" << endl;
    volScalarField Sal
    (
        IOobject
        (
            "Sal",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading sediment particle field Sed\n" << endl;
    volScalarField Sed
    (
        IOobject
        (
            "Sed",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Settling velocity field Vs\n" << endl;
    volVectorField Vs
    (
        IOobject
        (
            "Vs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0.0", U.dimensions(), vector(0.0, 0.0, 0.0)),
        U.boundaryField().types()
    );

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    dimensionedVector V0
    (
        transportProperties.lookup("V0")
    );

    Vs=V0;;

    #include "createPhi.H"

    #include "readTransportProperties.H"

    Info<< "Creating turbulence model\n" << endl;

    //New: general turbulence models (RAS, LES, Laminar)
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    //Original, only for RAS models
//    autoPtr<incompressible::RASModel> turbulence
//    (
//        incompressible::RASModel::New(U, phi, laminarTransport)
//    );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        calc_density(Sal,Sed)
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        //p_rgh = p - rhok*gh;
    }
