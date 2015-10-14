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

Class
    Foam::incompressible::RASModels::UQkOmega2006

Group
    grpIcoRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible flows.

    References:
	http://turbmodels.larc.nasa.gov/wilcox.html

	Turbulence Modeling for CFD (3rd Edition), David C. Wilcox, 2006

SourceFiles
    UQkOmega2006.C

\*---------------------------------------------------------------------------*/

#ifndef UQkOmega2006_H
#define UQkOmega2006_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class UQkOmega2006C2 Declaration
\*---------------------------------------------------------------------------*/

class UQkOmega2006
:
    public RASModel
{

protected:

    // Protected data

        // Model coefficients
			     IOdictionary modelCoefficients_;
           IOdictionary outputModelCoefDict_;
            dimensionedScalar Cmu_;
            //dimensionedScalar beta_; 				Commented for blending function
            dimensionedScalar beta0_;
            //dimensionedScalar kappa_;
            dimensionedScalar alpha_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;
	    	    dimensionedScalar Clim_;
	    	    dimensionedScalar alphadScalar_;

        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField nut_;
	    	volScalarField fBeta_;
	    	volScalarField Chi_;
	    	volScalarField absChi_;
	    	volScalarField beta_;
            volScalarField alphad_;

        // Protected Member Functions
        	dimensionedScalar alpha() const;


public:

    //- Runtime type information
    TypeName("UQkOmega2006");

    // Constructors

        //- Construct from components
        UQkOmega2006
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~UQkOmega2006()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                //new volScalarField("DkEff", alphaK_*nut_ + nu())
		new volScalarField("DkEff", alphaK_*k_/omega_ + nu())
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                //new volScalarField("DomegaEff", alphaOmega_*nut_ + nu())
		new volScalarField("DomegaEff", alphaOmega_*k_/omega_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;


        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
