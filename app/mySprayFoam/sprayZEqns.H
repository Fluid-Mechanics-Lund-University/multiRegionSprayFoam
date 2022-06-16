{
    tmp<fv::convectionScheme<scalar> > mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
	
	if(semiImplicit)
	{
		tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Z, dimMass/dimTime));
		fvScalarMatrix& fvm = tfvm.ref();
		forAll(Y,i)
		{
			fvm += parcels.SYi(i, Z);
		}
		fvScalarMatrix ZEqn
		(
			fvm::ddt(rho, Z)
		  + mvConvection->fvmDiv(phi, Z)
		  - fvm::laplacian((turbulence->mu()/Sc+turbulence->mut()/Sct), Z)
			==
			 tfvm
		  + fvOptions(rho,Z)
		);
		ZEqn.relax();
		fvOptions.constrain(ZEqn);
		ZEqn.solve("Yi");
		fvOptions.correct(Z);		
	}
	else
	{
		fvScalarMatrix ZEqn
		(
			fvm::ddt(rho, Z)
		  + mvConvection->fvmDiv(phi, Z)
		  - fvm::laplacian((turbulence->mu()/Sc+turbulence->mut()/Sct), Z)
			==
			 parcels.Srho(Z)
		  + fvOptions(rho,Z)
		);
		ZEqn.relax();
		fvOptions.constrain(ZEqn);
		ZEqn.solve("Yi");
		fvOptions.correct(Z);		
	}

	
    Chi = Cchi*turbulence->epsilon()/max(turbulence->k(),dimensionedScalar("smallK",dimVelocity*dimVelocity,SMALL))*Z2;
    Chi.correctBoundaryConditions();
	
    fvScalarMatrix Z2Eqn
    (
        fvm::ddt(rho, Z2)
      + mvConvection->fvmDiv(phi, Z2)
      - fvm::laplacian(turbulence->mu()/Sc+turbulence->mut()/Sct, Z2)
      ==
      + 2.0 * (turbulence->mu()/Sc+turbulence->mut()/Sct) *sqr(mag(fvc::grad(Z)))
      - rho*Chi
    );
	Z2Eqn.relax();
	fvOptions.constrain(Z2Eqn);
	Z2Eqn.solve("Yi");
	fvOptions.correct(Z2);

}
