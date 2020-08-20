/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    gravityDensityCurrentFoam based on buoyantBoussinesqPimpleFoam

    It accounts for density difference due to both scalar (such as
    salt and temperature) and particles (such as sediment; only one
    particle size is allowed right now).


Description
    Transient solver for buoyant, turbulent flow of incompressible fluids

    Uses the Boussinesq approximation.

Author
    Xiaofeng Liu, Penn State University
 
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "fvIOoptionList.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//constitutive relationship between excessive density and sal/sed f(Sal,Sed)
////rho/rho_0= F(Sal,Sed)
//Assume the units: Sal (PPT), Sed (kg/m^3)
tmp<volScalarField> calc_density(const volScalarField& Sal, const volScalarField& Sed)
{
   //temperature in degC
   scalar T = 20.0;

   //reference density 
   scalar rhow = 999.842594 + 6.793952e-2*T
                -9.095290e-3*pow(T,2)+1.001685e-4*pow(T,3)
                -1.120083e-6*pow(T,4)+6.536332e-9*pow(T,5);

   return 1.0+( (8.24493e-1-4.0899e-3*T+7.6438e-5*pow(T,2)-8.2467e-7*pow(T,3)+5.3875e-9*pow(T,4))*Sal
               +(-5.72466e-3+1.0227e-4*T-1.6546e-6*pow(T,2))*pow(Sal,1.5)
               +4.8314e-4*pow(Sal,2)
               + Sed
              )/rhow; 

   //return 1.0+1.65*Sed+0.001*Sal;
   //   return 1.0+0.001*S;
   //   return 1.0+0.000774*S;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "SEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
