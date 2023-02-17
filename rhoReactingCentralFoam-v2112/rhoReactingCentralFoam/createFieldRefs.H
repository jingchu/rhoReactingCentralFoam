volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}

PtrList<volScalarField> DiffsY(Y.size()); 
forAll(Y, i)
{
    DiffsY.set
    (
        i,
	new volScalarField 
	(
	    IOobject
    	    (
                "DiffsY",
                mesh.time().timeName(),
         	mesh,
         	IOobject::NO_READ,
         	IOobject::NO_WRITE
     	    ),
            mesh,
            dimensionedScalar("zero", mu.dimensions(), 0.0),
            zeroGradientFvPatchScalarField::typeName
	)
    );    
}
PtrList<volScalarField> HsY(Y.size()); 
forAll(Y, i)
{
    HsY.set
    (
        i,
	new volScalarField 
	(
	    IOobject
    	    (
                "HsY",
                mesh.time().timeName(),
         	mesh,
         	IOobject::NO_READ,
         	IOobject::NO_WRITE
     	    ),
            mesh,
            dimensionedScalar("zero", e.dimensions(), 0.0),
            zeroGradientFvPatchScalarField::typeName
	)
    );    
}

