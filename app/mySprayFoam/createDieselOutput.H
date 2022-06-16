fileName XYpath = "XY"+args.globalCaseName();
mkDir(XYpath);

OFstream Pav(XYpath/"p."+ runTime.timeName());
OFstream Tav(XYpath/"T."+ runTime.timeName());
OFstream kav(XYpath/"k."+ runTime.timeName());
OFstream epsav(XYpath/"eps."+ runTime.timeName());
OFstream sp95(XYpath/"sp95."+ runTime.timeName());
OFstream sp96(XYpath/"sp96."+ runTime.timeName());
OFstream sp97(XYpath/"sp97."+ runTime.timeName());
OFstream sp98(XYpath/"sp98."+ runTime.timeName());
OFstream sp99(XYpath/"sp99."+ runTime.timeName());
OFstream sp99_9(XYpath/"sp99.9."+ runTime.timeName());
OFstream spPVL0(XYpath/"spPVL0."+ runTime.timeName());
OFstream spPVL1(XYpath/"spPVL1."+ runTime.timeName());
OFstream spVap1(XYpath/"vapPen1."+ runTime.timeName());

volScalarField localLVF
(
   IOobject
   (
	  "localLVF",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
   ),
   mesh,
   dimensionedScalar("localLVF",dimensionSet (0,0,0,0,0,0,0),0.0)
);

// create projected cell list
const dictionary& subModelDict = parcels.subModelProperties();
const vector posInjector = subModelDict.subDict
(
    "injectionModels"
).subDict("model1").lookup("position");

vector dirInjection = subModelDict.subDict
(
    "injectionModels"
).subDict("model1").lookup("direction");
dirInjection /= mag(dirInjection);



