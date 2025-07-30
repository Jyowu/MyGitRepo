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
#line 31 "/home/jyothis/OneDrive/Fromserver/Hyrocyclone/Curvature_correction/Coarse/realKE/ClearCodes/realKE+kBuoyCurvS+eBuoy/constant/fvModels/kBouyancyCurvatureSource"
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
    // SHA1 = 21425ac7e2d8fb4d06d02127075c3f4f7d809b13
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void kBouyancyCurvatureSource_21425ac7e2d8fb4d06d02127075c3f4f7d809b13(bool load)
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

defineTypeNameAndDebug(kBouyancyCurvatureSourceFvModelscalarSource, 0);

addRemovableToRunTimeSelectionTable
(
    fvModel,
    kBouyancyCurvatureSourceFvModelscalarSource,
    dictionary
);


const char* const kBouyancyCurvatureSourceFvModelscalarSource::SHA1sum =
    "21425ac7e2d8fb4d06d02127075c3f4f7d809b13";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kBouyancyCurvatureSourceFvModelscalarSource::
kBouyancyCurvatureSourceFvModelscalarSource
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
        Info<<"construct kBouyancyCurvatureSource sha1: 21425ac7e2d8fb4d06d02127075c3f4f7d809b13"
            " from components\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

kBouyancyCurvatureSourceFvModelscalarSource::
~kBouyancyCurvatureSourceFvModelscalarSource()
{
    if (false)
    {
        Info<<"destroy kBouyancyCurvatureSource sha1: 21425ac7e2d8fb4d06d02127075c3f4f7d809b13\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kBouyancyCurvatureSourceFvModelscalarSource::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"kBouyancyCurvatureSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


void kBouyancyCurvatureSourceFvModelscalarSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"kBouyancyCurvatureSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 52 "/home/jyothis/OneDrive/Fromserver/Hyrocyclone/Curvature_correction/Coarse/realKE/ClearCodes/realKE+kBuoyCurvS+eBuoy/constant/fvModels/kBouyancyCurvatureSource"
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

            // Get the input data END



            // Get the required fields
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
               
		
            // Compute the source term
            const volScalarField::Internal Gb = nut/(ScT)*(g & fvc::grad(rho));

            const scalarField& V = mesh().V();

            Info<<endl<<endl<<"k Buoyancy Source computed="<<endl; 


  //#############################################BOUYANCY END ##############################################              



  //#############################################CURVATURE START ##############################################              



            tmp<volTensorField> tgradU = fvc::grad(U); // gradient of velocity field
            tmp<volSymmTensorField> tS = symm(tgradU()); // symmetric part of the gradient (strain rate tensor)
    
            tmp<volTensorField> tOmega = skew(tgradU()); // antisymmetric part of the gradient (vorticity tensor)
            volScalarField w= epsilon/(0.09* k + k_small); // turbulence dissipation rate divided by a small constant times k
   

            volScalarField S2("S2", 2.0 * magSqr(tS())); // square of the strain rate tensor
       
            volScalarField Omega2("Omega2", 2.0 * magSqr(tOmega()));
       
       
       
            volScalarField Omegamag = mag(tOmega()); // magnitude of the vorticity tensor

//^^^^^^^^^^^^^^^^^^calculating G^^^^^^^^^^^^^^^^^^^^^^


            const volScalarField::Internal G
            (
                "G",
                nut*(tgradU().v() && dev(twoSymm(tgradU().v())))
            );

            tgradU.clear();

//^^^^^^^^^^^^^^^^^^^^END calculating G^^^^^^^^^^^^^^^^^^^^^^


//^^^^^^^^^^^^^^^^^^^^^^^^^^^^calculating rStar^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            volScalarField rStar
            (
                "rStar",
                sqrt(S2) / ((sqrt(Omega2 ))+omega_small) // ratio of the square root of the strain rate tensor squared to the square root of the vorticity tensor squared, ensuring it does not go below a small constant
            );
// ratio of the square root of the strain rate tensor squared to the square root of the vorticity tensor squared, ensuring it does not go below a small constant

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^END calculating rStar^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^calculating rTilda^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

            volScalarField D2("D2", max(S2, 0.09 * sqr(w))); // maximum of the strain rate tensor squared and a small constant times the square of tefomega
            // This is a measure of the turbulence intensity, ensuring it does not go below a certain
            volScalarField Dtilda("Dtilda", pow(sqrt(D2), 3.0) * Omegamag); // cube of the square root of D2 plus a small constant, scaled by the magnitude of the vorticity tensor

           
            volSymmTensorField DSDt = fvc::ddt(tS())+fvc::div(phi, tS()); // time derivative of the strain rate tensor plus divergence of the flux of the strain rate tensor
            
            volTensorField DSDtfull("DSDtfull", (tOmega()&tOmega()) * 0.0 ); //Initialize DSDtfull with zeros


// Copy the components of DSDt to DSDtfull to ensure it DSDtfull is a full tensor field
            forAll(DSDtfull, celli)
            {
                    DSDtfull[celli].component(0) = DSDt[celli].component(0);
                    DSDtfull[celli].component(1) = DSDt[celli].component(1);
                    DSDtfull[celli].component(2) = DSDt[celli].component(2);
                    DSDtfull[celli].component(3) = DSDt[celli].component(1);
                    DSDtfull[celli].component(4) = DSDt[celli].component(3);
                    DSDtfull[celli].component(5) = DSDt[celli].component(4);
                    DSDtfull[celli].component(6) = DSDt[celli].component(2);
                    DSDtfull[celli].component(7) = DSDt[celli].component(4);
                    DSDtfull[celli].component(8) = DSDt[celli].component(5);
            }


            volTensorField rt1("rt1", 2.0 * (tOmega() & tS().T() ));// Calculate rTilda as a tensor field based on the vorticity tensor, strain rate tensor, and the production term
            volScalarField rTilda("rTilda", ((rt1 && DSDtfull) / ( Dtilda + dimensionedScalar("Dsmall", dimensionSet(0, 0, -4, 0, 0), SMALL))) ); // Calculate rTilda as a tensor field based on the vorticity tensor, strain rate tensor, and the production term

        
        
            tOmega.clear();
        
            tS.clear();



            scalar cr1 = 1.0;
            scalar cr2 = 2.0;
            scalar cr3 = 1.0;
            scalar C_curv = 1.0;


            volScalarField frotation =  ((scalar(1.0) + cr1)
                            *(scalar(2.0) * rStar / (scalar(1.0) + rStar+SMALL)))
                            * (scalar(1.0) - cr3 * atan(cr2 * rTilda+SMALL))
                            - cr1;       // Calculate the rotation function based on the strain rate and vorticity tensors


            volScalarField fr1 = max(min(frotation, scalar(1.25)), scalar(0.0)); // Ensure the rotation function is within the range [0, 1.25]




            volScalarField fr2 = max(scalar(0.0), scalar(1.0)+ C_curv*(fr1-scalar(1.0))); // Calculate the curvature correction factor based on the rotation function and a curvature constant



            scalarField& kSource = eqn.source();

            // Add the Buouyancy source term to the equation
            kSource += Gb*V;

            // Add the curvature correction source term to the equation
            kSource += (fr2 - 1.0) *rho()*G*V;
//}}} end code
}


void kBouyancyCurvatureSourceFvModelscalarSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (false)
    {
        Info<<"kBouyancyCurvatureSourceFvModelscalarSource::addSup()\n";
    }

//{{{ begin code
    
//}}} end code
}


bool kBouyancyCurvatureSourceFvModelscalarSource::movePoints()
{
    set_.movePoints();
    return true;
}


void kBouyancyCurvatureSourceFvModelscalarSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void kBouyancyCurvatureSourceFvModelscalarSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void kBouyancyCurvatureSourceFvModelscalarSource::distribute
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

