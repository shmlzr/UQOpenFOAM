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

\*---------------------------------------------------------------------------*/

#include "UQLaunderSharma.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(UQLaunderSharma, 0);
addToRunTimeSelectionTable(RASModel, UQLaunderSharma, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> UQLaunderSharma::fMu() const
{
    return exp(-3.4/sqr(scalar(1) + sqr(k_)/(nu()*epsilonTilda_)/50.0));
}


tmp<volScalarField> UQLaunderSharma::f2() const
{
    return
        scalar(1)
      - 0.3*exp(-min(sqr(sqr(k_)/(nu()*epsilonTilda_)), scalar(50.0)));
}

// Ceps1 and sigmaEps are chosen according to Edeling et al. (2015)
dimensionedScalar UQLaunderSharma::Ceps1() const
{
    return
	(Ceps2_-1.0)/2.09 + 1.0;

}
dimensionedScalar UQLaunderSharma::sigmaEps() const
{
    return
	pow(kappa_, 2.0)/(pow(Cmu_, 1.0/2.0)*(Ceps2_ - Ceps1()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

UQLaunderSharma::UQLaunderSharma
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    modelCoefficients_
    (IOdictionary
        (
        	IOobject
        	(
        	    "modelCoefficients",
        	    runTime_.constant(),
        	    mesh_,
        	    IOobject::MUST_READ_IF_MODIFIED,
        	    IOobject::NO_WRITE
        	)
    	  )
     ),
     outputModelCoefDict_
     (IOdictionary
       (
          IOobject
         	(
         	    "outputModelCoefDict",
         	    runTime_.constant(),
         	    mesh_,
         	    IOobject::READ_IF_PRESENT,
         	    IOobject::AUTO_WRITE
         	)
        )
    ),
    Cmu_
    (
	     readScalar ( modelCoefficients_.lookup("Cmu") )
    ),
    Ceps2_
    (
	     readScalar ( modelCoefficients_.lookup("Ceps2") )
    ),
    sigmaK_
    (
	     readScalar ( modelCoefficients_.lookup("Sigmak") )
    ),
    kappa_
    (
	     readScalar ( modelCoefficients_.lookup("Kappa") )
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    epsilonTilda_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateLowReNut("nut", mesh_)
    )

{
    bound(k_, kMin_);
    bound(epsilonTilda_, epsilonMin_);

    nut_ = Cmu_*fMu()*sqr(k_)/epsilonTilda_;
    nut_.correctBoundaryConditions();

    //printCoeffs();

    Info << "\n=================================================================================" << endl;
    Info <<   "==================== C L O S U R E   C O E F F I C I E N T S ====================" << endl;
    Info << "Ceps2    = " << Ceps2_.value() << endl;
    Info << "Cmu      = " << Cmu_.value() << endl;
    Info << "Sigmak   = " << sigmaK_.value() << endl;
    Info << "Kappa    = " << kappa_.value() << endl;
    Info << "Ceps1    = " << Ceps1().value() << endl;
    Info << "SigmaEps = " << sigmaEps().value() << endl;
    Info << "===================================================================================\n" << endl;

    outputModelCoefDict_.add< scalar >("Ceps2", Ceps2_.value(),true);
    outputModelCoefDict_.add< scalar >("Cmu", Cmu_.value(),true);
    outputModelCoefDict_.add< scalar >("Sigmak", sigmaK_.value(),true);
    outputModelCoefDict_.add< scalar >("Kappa", kappa_.value(),true);
    outputModelCoefDict_.add< scalar >("Ceps1", Ceps1().value(),true);
    outputModelCoefDict_.add< scalar >("SigmaEps", sigmaEps().value(),true);
    outputModelCoefDict_.regIOobject::write();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> UQLaunderSharma::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> UQLaunderSharma::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> UQLaunderSharma::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> UQLaunderSharma::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool UQLaunderSharma::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void UQLaunderSharma::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volScalarField> S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G(GName(), nut_*S2);

    const volScalarField E(2.0*nu()*nut_*fvc::magSqrGradGrad(U_));
    const volScalarField D(2.0*nu()*magSqr(fvc::grad(sqrt(k_))));


    Info<< nl << "Cmu__" << Cmu_.value() << ", "
    	<< "Ceps2__" << Ceps2_.value() << ", "
	<< "sigmaK__" << sigmaK_.value() << ", "
    	<< "kappa__" << kappa_.value() << nl << endl;

    // Dissipation rate equation

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilonTilda_)
      + fvm::div(phi_, epsilonTilda_)
      - fvm::laplacian(DepsilonEff(), epsilonTilda_)
     ==
        Ceps1()*G*epsilonTilda_/k_
      - fvm::Sp(Ceps2_*f2()*epsilonTilda_/k_, epsilonTilda_)
      + E
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilonTilda_, epsilonMin_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp((epsilonTilda_ + D)/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ == Cmu_*fMu()*sqr(k_)/epsilonTilda_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
