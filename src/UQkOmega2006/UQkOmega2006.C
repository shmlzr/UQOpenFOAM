/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "UQkOmega2006.H"
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

defineTypeNameAndDebug(UQkOmega2006, 0);
addToRunTimeSelectionTable(RASModel, UQkOmega2006, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

// alpha are chosen according to Edeling et al. (2015)
/*dimensionedScalar UQkOmega2006::alpha() const
{
	return
	(beta0_/Cmu_) - (pow(kappa_, 2.0)/ (2*pow(Cmu_, 1.0/2.0)));
}*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

UQkOmega2006::UQkOmega2006
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
    	readScalar ( modelCoefficients_.lookup("Betas") )
    ),
    beta0_
    (
	    readScalar ( modelCoefficients_.lookup("Beta") )
    ),
		alpha_
		(
    	readScalar ( modelCoefficients_.lookup("Alpha") )
    ),
    alphaK_
    (
    	readScalar ( modelCoefficients_.lookup("Sigmas") )
    ),
    alphaOmega_
    (
	    readScalar ( modelCoefficients_.lookup("Sigma") )
    ),
    Clim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clim",
            coeffDict_,
            0.875 //Clim=7/8
        )
    ),
    alphadScalar_
    (
			dimensioned<scalar>::lookupOrAddToDict
			(
					"Sigmado",
					coeffDict_,
					0.125
			)
    ),
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
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
        autoCreateNut("nut", mesh_)
    ),
   fBeta_
    (
        IOobject
        (
            "fBeta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
   Chi_
    (
        IOobject
        (
            "Chi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_, dimless
    ),
   absChi_
    (
        IOobject
        (
            "absChi",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_, dimless
    ),
   beta_
    (
        IOobject
        (
            "beta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimless
    ),
   alphad_
    (
        IOobject
        (
            "alphad",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_, dimensionedScalar("zero", dimless, alphadScalar_.value())
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

   //nut_ = k_/omega_; //Standard OpenFOAM definition
    nut_ = k_/ max(omega_, Clim_*sqrt(2.0/Cmu_*magSqr(symm(fvc::grad(U_)))));; //FIXME: Replace 0.09 by betastar
    nut_.correctBoundaryConditions();


    //printCoeffs();

		Info << "\n=================================================================================" << endl;
		Info <<   "==================== C L O S U R E   C O E F F I C I E N T S ====================" << endl;
		Info << "Alpha   = " << alpha_.value() << endl;
		Info << "Beta    = " << beta0_.value() << endl;
		Info << "Betas   = " << Cmu_.value() << endl;
		Info << "Sigma   = " << alphaOmega_.value() << endl;
		Info << "Sigmas  = " << alphaK_.value() << endl;
		Info << "Clim    = " << Clim_.value() << endl;
		Info << "Sigmado = " << alphadScalar_.value() << endl;
		Info << "===================================================================================\n" << endl;

    outputModelCoefDict_.add< scalar >("Alpha", alpha_.value(),true);
    outputModelCoefDict_.add< scalar >("Beta", beta0_.value(),true);
    outputModelCoefDict_.add< scalar >("Betas", Cmu_.value(),true);
    outputModelCoefDict_.add< scalar >("Sigma", alphaOmega_.value(),true);
    outputModelCoefDict_.add< scalar >("Sigmas", alphaK_.value(),true);
    outputModelCoefDict_.add< scalar >("Clim", Clim_.value(),true);
    outputModelCoefDict_.add< scalar >("Sigmado", alphadScalar_.value(),true);
    outputModelCoefDict_.regIOobject::write();

    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> UQkOmega2006::R() const
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


tmp<volSymmTensorField> UQkOmega2006::devReff() const
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


tmp<fvVectorMatrix> UQkOmega2006::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> UQkOmega2006::divDevRhoReff
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


bool UQkOmega2006::read()
{
    if (RASModel::read())
    {
    	//beta0_.readIfPresent(coeffDict());			//Beta
        //Cmu_.readIfPresent(coeffDict()); 			//Betas
        //beta_.readIfPresent(coeffDict());  		Must be commented for blending function
        //alphaK_.readIfPresent(coeffDict());			//Sigmas
        //alphaOmega_.readIfPresent(coeffDict());		//Sigma
		//alphaD_.readIfPresent(coeffDict());			//Sigmado
		//kappa_.readIfPresent(coeffDict());			//Kappa


        return true;
    }
    else
    {
        return false;
    }
}


void UQkOmega2006::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }



    volTensorField GradU(fvc::grad(U_));
    volSymmTensorField Sij(symm(GradU));
    volTensorField Omij(-skew(GradU));
    volScalarField StressLim(Clim_*sqrt(2.0/Cmu_)*mag(Sij));
    volSymmTensorField tauij(2.0*nut_*Sij-((2.0/3.0)*I)*k_);
    volVectorField Gradk(fvc::grad(k_));
    volVectorField Gradomega(fvc::grad(omega_));

    //volScalarField G(type() + ".G", tauij && GradU);
    volScalarField G(GName(), tauij && GradU); //for newer OF versions


    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

//START NEW STUFF FOR 2006.....................................

volScalarField alphadCheck_(Gradk & Gradomega);  //condition to change alphad_

forAll(alphad_,celli)
{
	if (alphadCheck_[celli] <= 0.0001)
	{
	alphad_[celli]=scalar(0);
	}else
	{
	alphad_[celli]=alphadScalar_.value(); //FIXME: sigmaDo
	}
}


volScalarField CDkOmega(alphad_/omega_*(Gradk & Gradomega)); //last term in NASA equations


Chi_ = (Omij & Omij) && Sij /pow((Cmu_*omega_),3);
absChi_ = mag(Chi_);
//fBeta_ = 1.0; //2D version. This term should be (1.0+85.0*absChi_)/(1.0+100.0*absChi_); for 3D
fBeta_ = (1.0+85.0*absChi_)/(1.0+100.0*absChi_); //3D
beta_ = beta0_*fBeta_; //FIXME: 0.0708 --> beta0



    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
        alpha_*G*omega_/k_
      - fvm::Sp(beta_*omega_, omega_)
      + CDkOmega //Crossflow diffusion term to match 2006
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(Cmu_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    //nut_ = k_/omega_; //Standard OpenFOAM definition
    nut_ = k_/ max(omega_, Clim_*sqrt(2.0/Cmu_*magSqr(symm(fvc::grad(U_)))));;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
