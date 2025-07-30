/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "codedFvModelTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "read.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
#line 38 "/home/jyothis/OneDrive/Fromserver/Hyrocyclone/Curvature_correction/Coarse/C_bouyantreal_ke_bouyant/constant/fvModels/kBouyancySource"
#include "fvm.H"
//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 1ce7d5adc0a79fa242b8a1050f1bb70fe03021f4
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void kBouyancySource_1ce7d5adc0a79fa242b8a1050f1bb70fe03021f4(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kBouyancySourceFvModelscalarSource, 0);

addRemovableToRunTimeSelectionTable
(
    fvModel,
    kBouyancySourceFvModelscalarSource,
    dictionary
);


const char* const kBouyancySourceFvModelscalarSource::SHA1sum =
    "1ce7d5adc0a79fa242b8a1050f1bb70fe03021f4";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kBouyancySourceFvModelscalarSource::
kBouyancySourceFvModelscalarSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    set_(mesh, coeffs())
{
    if (false)
    {
        Info<<"construct kBouyancySource sha1: 1ce7d5adc0a79fa242b8a1050f1bb70fe03021f4"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kBouyancySourceFvModelscalarSource::
~kBouyancySourceFvModelscalarSource()
{
    if (false)
    {
        Info<<"destroy kBouyancySource sha1: 1ce7d5adc0a79fa242b8a1050f1bb70fe03021f4\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kBouyancySourceFvModelscalarSource::addSup
(
    const VolField<scalar>& field,
    fvMatrix<scalar>& eqn
) const
{
    if (false)
    {
        Info<<"kBouyancySourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 48 "/home/jyothis/OneDrive/Fromserver/Hyrocyclone/Curvature_correction/Coarse/C_bouyantreal_ke_bouyant/constant/fvModels/kBouyancySource"
Info << "kCurvatureCorrection model constructed and activated." << endl;
             // Get the input data
             dimensionedVector g
             (
                 dimensionedVector(
                     "g",
                     dimAcceleration,
		     coeffs().lookup("g")
                 )
             ); 
             scalar ScT(readScalar(coeffs().lookup("ScT")));

             // Get the required fields
             const volScalarField& density = mesh().lookupObject<volScalarField>("rho");
             const volScalarField& nut = mesh().lookupObject<volScalarField>("nut");

	     // Get cell volume
	     const scalarField& V = mesh().V();
		
             // Compute the source term
             const scalarField Gb = nut.internalField()/ScT*(g & fvc::grad(density))*V;

             // Add the source term to the equation
             scalarField& kSource = eqn.source();
             kSource += Gb;
//}}} end code
}


void kBouyancySourceFvModelscalarSource::addSup
(
    const volScalarField& rho,
    const VolField<scalar>& field,
    fvMatrix<scalar>& eqn
) const
{
    if (false)
    {
        Info<<"kBouyancySourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


void kBouyancySourceFvModelscalarSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<scalar>& field,
    fvMatrix<scalar>& eqn
) const
{
    if (false)
    {
        Info<<"kBouyancySourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


bool kBouyancySourceFvModelscalarSource::movePoints()
{
    set_.movePoints();
    return true;
}


void kBouyancySourceFvModelscalarSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void kBouyancySourceFvModelscalarSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void kBouyancySourceFvModelscalarSource::distribute
(
    const polyDistributionMap& map
)
{
    set_.distribute(map);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace fv

// ************************************************************************* //

