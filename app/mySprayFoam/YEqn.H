tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    combustion->correct();
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
	   volScalarField& Yi = Y[i];

		if (i != inertIndex && composition.active(i))
		{
			fvScalarMatrix YEqn
			(
				fvm::ddt(rho, Yi)
			  + mvConvection->fvmDiv(phi, Yi)
			  - fvm::laplacian(turbulence->mut()/Sct + turbulence->mu()/Sc, Yi)
			 ==
				parcels.SYi(i, Yi)
			  + combustion->R(Yi)
			  + fvOptions(rho, Yi)
			);

			YEqn.relax();

			fvOptions.constrain(YEqn);

			YEqn.solve("Yi");

			fvOptions.correct(Yi);

			Yi.max(0.0);
			Yt += Yi;
		}									
	}
	Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}

