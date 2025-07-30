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
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude
#line 245 "/home/jyothis/OneDrive/Fromserver/Hyrocyclone/Curvature_correction/Coarse/realKE/ClearCodes/realKE+kBuoyCurvS+eBuoy/constant/fvModels/epsilonBouyancySource"
#include "fvm.H"
            #include "fvc.H"
            #include <cmath>
            #include "uniformDimensionedFields.H"
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
    // SHA1 = f5a038ec757c1512d68b57b1a19e6d7ea71d3fcd
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void epsilonBouyancySource_f5a038ec757c1512d68b57b1a19e6d7ea71d3fcd(bool load)
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

defineTypeNameAndDebug(epsilonBouyancySourceFvModelscalarSource, 0);

addRemovableToRunTimeSelectionTable
(
    fvModel,
    epsilonBouyancySourceFvModelscalarSource,
    dictionary
);


const char* const epsilonBouyancySourceFvModelscalarSource::SHA1sum =
    "f5a038ec757c1512d68b57b1a19e6d7ea71d3fcd";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonBouyancySourceFvModelscalarSource::
epsilonBouyancySourceFvModelscalarSource
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
        Info<<"construct epsilonBouyancySource sha1: f5a038ec757c1512d68b57b1a19e6d7ea71d3fcd"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

epsilonBouyancySourceFvModelscalarSource::
~epsilonBouyancySourceFvModelscalarSource()
{
    if (false)
    {
        Info<<"destroy epsilonBouyancySource sha1: f5a038ec757c1512d68b57b1a19e6d7ea71d3fcd\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonBouyancySourceFvModelscalarSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"epsilonBouyancySourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


void epsilonBouyancySourceFvModelscalarSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"epsilonBouyancySourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 258 "/home/jyothis/OneDrive/Fromserver/Hyrocyclone/Curvature_correction/Coarse/realKE/ClearCodes/realKE+kBuoyCurvS+eBuoy/constant/fvModels/epsilonBouyancySource"
//#############################################BOUYANCY START ##############################################  
           

            // Get the input data START

            const uniformDimensionedVectorField& g = mesh().lookupObject<uniformDimensionedVectorField>("g");

            IOdictionary momentumTransportDict
            (
                IOobject
                (
                    "momentumTransport",
                    mesh().time().constant(),
                    mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            scalar ScT = readScalar(momentumTransportDict.lookup("ScT"));

            const volVectorField& U = mesh().lookupObject<volVectorField>("U");
            const volScalarField& nut = mesh().lookupObject<volScalarField>("nut");
            const volScalarField& k    = mesh().lookupObject<volScalarField>("k");
            const volScalarField& epsilon = mesh().lookupObject<volScalarField>("epsilon");
            const surfaceScalarField& phi = mesh().lookupObject<surfaceScalarField>("phi");
            const volScalarField& alphas = mesh().lookupObject<volScalarField>("alpha.sludge");
            const volScalarField& alphaw = mesh().lookupObject<volScalarField>("alpha.water");


            // Create dimensioned SMALL constants
            dimensionedScalar rho_small("rho_small", dimMass/dimVolume, SMALL);
            dimensionedScalar k_small("k_small", dimEnergy/dimMass, SMALL);
            dimensionedScalar omega_small("omega_small", dimless/dimTime, SMALL);
               
            // Get cell volume
	        const scalarField& V = mesh().V();
            scalar C1 = 1.44;

            vector gHat(g.value()/mag(g.value()));
            volScalarField v(gHat & U);
            volScalarField u
            (
                mag(U- gHat*v)
                + dimensionedScalar(dimVelocity, small)
            );

             
		
            // Compute the source term
             const volScalarField::Internal Gb = nut/(ScT)*(g & fvc::grad(rho));


             // Compute the buoyancy term
        
             // Add the source term to the equation
             scalarField& epsilonSource = eqn.source();
             epsilonSource += C1*epsilon/(k + k_small)*tanh(mag(v)/u)*Gb*V;

             //Info<<"Source added ="<<endl<<endl;
//}}} end code
}


void epsilonBouyancySourceFvModelscalarSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"epsilonBouyancySourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


bool epsilonBouyancySourceFvModelscalarSource::movePoints()
{
    set_.movePoints();
    return true;
}


void epsilonBouyancySourceFvModelscalarSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void epsilonBouyancySourceFvModelscalarSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void epsilonBouyancySourceFvModelscalarSource::distribute
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

