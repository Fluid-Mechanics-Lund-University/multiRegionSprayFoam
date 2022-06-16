Info << "max T = " << max(T).value() << " K" << endl;
Info << "min T = " << min(T).value() << " K" << endl;
Info << "max p = " << max(p).value()/1.0e5 << " bar" << endl;
Info << "min p = " << min(p).value()/1.0e5 << " bar" << endl;
Info << "max rho = " << max(rho).value() << " kg/m3" << endl;
Info << "min rho = " << min(rho).value() << " kg/m3" << endl;

Pav << runTime.timeName() << "   " << p.weightedAverage(mesh.V()).value() << endl;
Tav << runTime.timeName() << "   " << T.weightedAverage(mesh.V()).value() <<"   " << max(T).value() << "   " << min(T).value()<<endl;

kav << runTime.timeName() << "   " << (turbulence->k())().weightedAverage(mesh.V()).value() << endl;
epsav << runTime.timeName() << "   " << (turbulence->epsilon())().weightedAverage(mesh.V()).value() << endl;
sp95 <<runTime.timeName() << " " << parcels.penetration(0.95) <<endl;
sp96 <<runTime.timeName() << " " << parcels.penetration(0.96) <<endl;
sp97 <<runTime.timeName() << " " << parcels.penetration(0.97) <<endl;
sp98 <<runTime.timeName() << " " << parcels.penetration(0.98) <<endl;
sp99 <<runTime.timeName() << " " << parcels.penetration(0.99) <<endl;
sp99_9 <<runTime.timeName() << " " << parcels.penetration(0.999) <<endl;

//https://ecn.sandia.gov/diesel-spray-combustion/computational-method/modeling-standards/
// Liquid penetration, projected liquid volume 0.2e-3 mm3/mm2, 2e-3 mm3/mm2
//
scalar maxVapPen1 = -GREAT;

forAll(T, celli)
{
	vector posCells = (T.mesh().cellCentres()[celli] - posInjector);
//	if(mag(posCells/mag(posCells) & dirInjection/mag(dirInjection)) > 0.98)

	if(mag(posCells/mag(posCells) & dirInjection/mag(dirInjection)) > 0.9)
	{	
		if((mag(posCells&dirInjection) > maxVapPen1) && (Z[celli] > 0.001)) 
		{
			maxVapPen1 = mag(posCells&dirInjection);
		}
	}
}

reduce(maxVapPen1, maxOp<scalar>());
Info <<" maxVapPen1 "<< " "<< maxVapPen1 << endl;
if(maxVapPen1 < 0)
{
	maxVapPen1 = 0.0;
}

spVap1 << runTime.timeName() << " " << maxVapPen1 << endl;



