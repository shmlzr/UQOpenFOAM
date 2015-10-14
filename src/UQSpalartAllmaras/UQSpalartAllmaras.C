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

#include "UQSpalartAllmaras.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(UQSpalartAllmaras, 0);
addToRunTimeSelectionTable(RASModel, UQSpalartAllmaras, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> UQSpalartAllmaras::chi() const
{
    return nuTilda_/nu();
}


tmp<volScalarField> UQSpalartAllmaras::fv1(const volScalarField& chi) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> UQSpalartAllmaras::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
    /*if (ashfordCorrection_)
    {
        return 1.0/pow3(scalar(1) + chi/Cv2_);
    }
    else
    {
        return 1.0 - chi/(1.0 + chi*fv1);
    }*/
}


tmp<volScalarField> UQSpalartAllmaras::fv3
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "fv3",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("fv3", dimless, 1),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    /*
    if (ashfordCorrection_)
    {
        //const volScalarField chiByCv2((1/Cv2_)*chi);

        return
            (scalar(1) + chi*fv1)
           *(1/Cv2_)
           *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
           /pow3(scalar(1) + chiByCv2);
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fv3",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("fv3", dimless, 1),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }*/
}


tmp<volScalarField> UQSpalartAllmaras::fw(const volScalarField& Stilda) const
{
    volScalarField r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*d_)
            ),
            scalar(10.0)
        )
    );
    Info << " rMax = " << max(r.internalField()) << endl;

    r.boundaryField() == 0.0;

    const volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

UQSpalartAllmaras::UQSpalartAllmaras
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

    sigmaNut_
    (
        readScalar ( modelCoefficients_.lookup("Sigma") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )*/
    ),
    kappa_
    (
        readScalar ( modelCoefficients_.lookup("Kappa") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )*/
    ),

    Cb1_
    (
        readScalar ( modelCoefficients_.lookup("Cb1") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )*/
    ),
    Cb2_
    (
        readScalar ( modelCoefficients_.lookup("Cb2") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )*/
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        readScalar ( modelCoefficients_.lookup("Cw2") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )*/
    ),
    Cw3_
    (
        readScalar ( modelCoefficients_.lookup("Cw3") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )*/
    ),
    Cv1_
    (
        readScalar ( modelCoefficients_.lookup("Cv1") )
        /*dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )*/
    ),
    /*Cv2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv2",
            coeffDict_,
            5.0
        )
    ),*/

    ashfordCorrection_(coeffDict_.lookupOrDefault("ashfordCorrection", false)),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    d_(mesh_)
{
    //printCoeffs();

    Info << "\n=================================================================================" << endl;
    Info <<   "==================== C L O S U R E   C O E F F I C I E N T S ====================" << endl;
    Info << "Sigma   = " << sigmaNut_.value() << endl;
    Info << "Kappa   = " << kappa_.value() << endl;
    Info << "Cb1     = " << Cb1_.value() << endl;
    Info << "Cb2     = " << Cb2_.value() << endl;
    Info << "Cw1     = " << Cw1_.value() << endl;
    Info << "Cw2     = " << Cw2_.value() << endl;
    Info << "Cw3     = " << Cw3_.value() << endl;
    Info << "Cv1     = " << Cv1_.value() << endl;
    Info << "===================================================================================\n" << endl;

    outputModelCoefDict_.add< scalar >("Sigma", sigmaNut_.value(),true);
    outputModelCoefDict_.add< scalar >("Kappa", kappa_.value(),true);
    outputModelCoefDict_.add< scalar >("Cb1", Cb1_.value(),true);
    outputModelCoefDict_.add< scalar >("Cb2", Cb2_.value(),true);
    outputModelCoefDict_.add< scalar >("Cw1", Cw1_.value(),true);
    outputModelCoefDict_.add< scalar >("Cw2", Cw2_.value(),true);
    outputModelCoefDict_.add< scalar >("Cw3", Cw3_.value(),true);
    outputModelCoefDict_.add< scalar >("Cv1", Cv1_.value(),true);
    outputModelCoefDict_.regIOobject::write();


    if (ashfordCorrection_)
    {
        Info<< "    Employing Ashford correction" << endl;
    }
    else
    {
        Info<< "No Ashford correction, standard version of Spalart-Allmaras is used" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> UQSpalartAllmaras::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", (nuTilda_ + nu())/sigmaNut_)
    );
}


tmp<volScalarField> UQSpalartAllmaras::k() const
{
    WarningIn("tmp<volScalarField> UQSpalartAllmaras::k() const")
        << "Turbulence kinetic energy not defined for Spalart-Allmaras model. "
        << "Returning zero field" << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
        )
    );
}


tmp<volScalarField> UQSpalartAllmaras::epsilon() const
{
    WarningIn("tmp<volScalarField> UQSpalartAllmaras::epsilon() const")
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
        )
    );
}


tmp<volSymmTensorField> UQSpalartAllmaras::R() const
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
            ((2.0/3.0)*I)*k() - nut()*twoSymm(fvc::grad(U_))
        )
    );
}


tmp<volSymmTensorField> UQSpalartAllmaras::devReff() const
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


tmp<fvVectorMatrix> UQSpalartAllmaras::divDevReff(volVectorField& U) const
{
    const volScalarField nuEff_(nuEff());

    return
    (
      - fvm::laplacian(nuEff_, U)
      - fvc::div(nuEff_*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> UQSpalartAllmaras::divDevRhoReff
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


bool UQSpalartAllmaras::read()
{
    if (RASModel::read())
    {
        sigmaNut_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());

        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        //Cv2_.readIfPresent(coeffDict());

        ashfordCorrection_.readIfPresent("ashfordCorrection", coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void UQSpalartAllmaras::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        // Re-calculate viscosity
        nut_ = nuTilda_*fv1(this->chi());
        nut_.correctBoundaryConditions();

        return;
    }

    if (mesh_.changing())
    {
        d_.correct();
    }

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    const volScalarField Stilda
    (
        fv3(chi, fv1)*::sqrt(2.0)*mag(skew(fvc::grad(U_)))
      + fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_)
    );

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - fvm::laplacian(DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*magSqr(fvc::grad(nuTilda_))
     ==
        Cb1_*Stilda*nuTilda_
      - fvm::Sp(Cw1_*fw(Stilda)*nuTilda_/sqr(d_), nuTilda_)
    );

    nuTildaEqn().relax();
    solve(nuTildaEqn);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    // Re-calculate viscosity
    nut_.internalField() = fv1*nuTilda_.internalField();
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
