/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    rhoCentralFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on
    central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

    Chemistry reaction solver for multi-gas.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes."
        "Chemistry reaction solver for multi-gas."
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "createFvOptions.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (!LTS)
        {
            #include "setDeltaT.H"

            ++runTime;

            // Do any mesh changes
            mesh.update();
        }

        // --- Directed interpolation of primitive fields onto faces

        surfaceScalarField rho_pos(interpolate(rho, pos));
        surfaceScalarField rho_neg(interpolate(rho, neg));

        surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
        surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

        volScalarField rPsi("rPsi", 1.0/psi);
        surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
        surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));

        surfaceScalarField e_pos(interpolate(e, pos, T.name()));
        surfaceScalarField e_neg(interpolate(e, neg, T.name()));

        surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
        surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

        surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
        surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

	    PtrList<surfaceScalarField> Y_pos(Y.size()), Y_neg(Y.size()); 
        forAll(Y,i)
	    {
	        Y_pos.set(i,new surfaceScalarField (fvc::interpolate(Y[i], pos, "reconstruct(Y)")));
	        Y_neg.set(i,new surfaceScalarField (fvc::interpolate(Y[i], neg, "reconstruct(Y)")));
	    }

        surfaceScalarField phiv_pos("phiv_pos", U_pos & mesh.Sf());
        // Note: extracted out the orientation so becomes unoriented
        phiv_pos.setOriented(false);
        surfaceScalarField phiv_neg("phiv_neg", U_neg & mesh.Sf());
        phiv_neg.setOriented(false);

        // Make fluxes relative to mesh-motion
        if (mesh.moving())
        {
            surfaceScalarField meshPhi(mesh.phi());
            meshPhi.setOriented(false);
            phiv_pos -= meshPhi;
            phiv_neg -= meshPhi;
        }

        volScalarField c("c", sqrt(thermo.Cp()/thermo.Cv()*rPsi));
        surfaceScalarField cSf_pos
        (
            "cSf_pos",
            interpolate(c, pos, T.name())*mesh.magSf()
        );

        surfaceScalarField cSf_neg
        (
            "cSf_neg",
            interpolate(c, neg, T.name())*mesh.magSf()
        );

        surfaceScalarField ap
        (
            "ap",
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );

        surfaceScalarField am
        (
            "am",
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos("a_pos", ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf("aSf", am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
        surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "centralCourantNo.H"

        if (LTS)
        {
            #include "setRDeltaT.H"

            ++runTime;
        }

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiU(aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg);
        // Note: reassembled orientation from the pos and neg parts so becomes
        // oriented
        phiU.setOriented(true);

        surfaceVectorField phiUp(phiU + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

        surfaceScalarField phiEp
        (
            "phiEp",
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );

        // Make flux for pressure-work absolute
        if (mesh.moving())
        {
            surfaceScalarField meshPhi(mesh.phi());
            meshPhi.setOriented(false);
            phiEp += meshPhi*(a_pos*p_pos + a_neg*p_neg);
        }

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() =
            rhoU()
           /rho();
        fvOptions.correct(U);			
	U.max(U_min);
	U.min(U_max);
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

		Info<< "min/max(U) = " << min(mag(U)).value() << ", " << max(mag(U)).value() << endl;
		
        // --- Solve species equation        
		PtrList<surfaceScalarField> phiY(Y.size()); 
	    forAll(Y, i) 
	    { 
	        phiY.set
	       (
		      i,
		      new surfaceScalarField
		      ( 
		         (aphiv_pos*(rho_pos*Y_pos[i]) + aphiv_neg*(rho_neg*Y_neg[i]))	 
		      )	
	        ); 
	    }
		
	    Yt=0.0;  
	    forAll(Y, i)
	    {
			if (Y[i].name() != inertSpecie)
		    {   
				solve(fvm::ddt(rhoY[i]) + fvc::div(phiY[i]));
				Y[i] = rhoY[i]/rho;  
				Y[i].correctBoundaryConditions();
				Yt += Y[i];
				Y[i].min(1.0);
		    }
		}
		Yt.correctBoundaryConditions();
         	Info<< "min/max(Yt after convection) = " << min(Yt).value() << ", " << max(Yt).value() << endl;
	
		Y[inertIndex] =1.0- Yt;
		Y[inertIndex].correctBoundaryConditions();
        
		forAll(Y, i)
	    {
            Y[i].correctBoundaryConditions();
	        rhoY[i].boundaryFieldRef() == rho.boundaryField()*Y[i].boundaryField();
	    }
	    if(!inviscid)
	    {
	        forAll(Y,i)
			{
			DiffsY[i]=turbulence->mut()/Sct+turbulence->mu()/Sc[i];
			}
			Y[inertIndex]=1.0;
	        forAll(Y, i)
	        {
		       if (Y[i].name() != inertSpecie)
		       {   
		          volScalarField& Yi = Y[i];
		          fvScalarMatrix YiEqn
		          (
			         fvm::ddt(rho, Yi)- fvc::ddt(rho, Yi) 
		             == fvm::laplacian(DiffsY[i], Yi)
		          );
		          YiEqn.relax();   
				  YiEqn.solve(mesh.solver("Yi"));
				  //Yi.max(0.0); 
                  Yi.min(1.0);    
				  Y[inertIndex] -=Y[i];
		       }
	        }
            Yt=0.0;    
			Yt.correctBoundaryConditions();
	        forAll(Y,i)	    
            {   
                Yt +=Y[i];
                Y[i].correctBoundaryConditions();  
                rhoY[i] = rho*Y[i];
            }
            Yt.correctBoundaryConditions();
	    }
        Info<< "min/max(Yt after diffusion) = " << min(Yt).value() << ", " << max(Yt).value() << endl;

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            "sigmaDotU",
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + fvc::dotInterpolate(mesh.Sf(), tauMC)
            )
          & (a_pos*U_pos + a_neg*U_neg)
        );

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );

        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() ==
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        if (!inviscid)
        {
			//sensible enthalpy diffusion term induced by species diffusion
			forAll(HsY,i)
			{
				forAll(HsY[i],celli)
				{
				   HsY[i][celli]=composition.Hs(i,p[celli],T[celli]);
				}
				HsY[i].correctBoundaryConditions();
			}
			volVectorField HsDiff("HsDiff",DiffsY[0]*HsY[0]*fvc::grad(Y[0]));
			
			for(label i=1;i<Y.size();i++)
			{
				HsDiff +=DiffsY[i]*HsY[i]*fvc::grad(Y[i]);	 
			}
			surfaceScalarField HsDiffTerm=(mesh.Sf()&(fvc::interpolate(HsDiff)))(); 
			
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvc::div(HsDiffTerm)
              - fvm::laplacian(turbulence->alphaEff(), e)
            );
	    fvOptions.correct(e);
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }
        Info<< "min/max(T) after convection & diffusion = " << min(T).value() << ", " << max(T).value() << endl;
        Info<< "min/max(p) after convection & diffusion = " << min(p).value() << ", " << max(p).value() << endl;

        if(SolveChemistry)
        {
            combustion->correct();
			Qdot = combustion->Qdot();	
			Qdot.correctBoundaryConditions();	
            Info << "Qdot = " << gSum(Qdot.internalField())<<endl;
			
            Yt = 0.0;   
            forAll(Y,i)
            {
				if(i!=inertIndex && composition.active(i))
				{
					solve(fvm::ddt(rhoY[i])-fvc::ddt(rhoY[i]) == combustion->R2(Y[i]));
                    Y[i]=rhoY[i]/rho;   
					Y[i].max(0.0);
				    Yt += Y[i];
				}
            }
            Yt.correctBoundaryConditions();
			Y[inertIndex] = scalar(1) - Yt;
			Y[inertIndex].max(0.0);
			
			Yt = 0.0;
			forAll(Y, i)
			{
				Yt += Y[i];
			}
            Info<<"min/max(Yt after reaction) = "<<min(Yt).value()<<" , "<<max(Yt).value()<<endl;
            forAll(Y, i)
	        {
            	Y[i].correctBoundaryConditions();
                rhoY[i].boundaryFieldRef() == rho.boundaryField()*Y[i].boundaryField();
            }
	
            /*fvScalarMatrix EEqn
	        (
                fvm::ddt(rhoE)-fvc::ddt(rhoE) == combustion->Qdot()
            );
			Info<<"breakpoint---relax---1---"<<endl;
           	EEqn.relax(); 
			Info<<"breakpoint---relax---2---"<<endl;
			EEqn.solve();*/
			//Info<<"breakpoint---relax---2---"<<endl;
			solve(fvm::ddt(rhoE)-fvc::ddt(rhoE) == combustion->Qdot());
			//Info<<"breakpoint---relax---3---"<<endl;
            e = rhoE/rho - 0.5*magSqr(U); 
            e.correctBoundaryConditions();
	    fvOptions.correct(e);
            thermo.correct();
            rhoE.boundaryFieldRef() == 
				rho.boundaryField()*
				(
					e.boundaryField() + 0.5*magSqr(U.boundaryField())
				);
        }
        p.ref() =
            rho()
           /psi();
	fvOptions.correct(p);
        p.correctBoundaryConditions();
        rho == psi*p;
        if(1)
	{							
	   rhoU = rho*U;
	   forAll(Y,i){rhoY[i] = rho*Y[i];}	
           rhoE=rho*(e+0.5*magSqr(U));	
	}

        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
