Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

Switch  semiImplicit(chemistryProperties.lookupOrDefault("semiImplicit",true));
label invStep(chemistryProperties.lookupOrDefault("invStep",10));		

// Prandlt : 0.7 - 0.9 (0.85)
dimensionedScalar Prt
(
   chemistryProperties.lookup("Prt")
);

dimensionedScalar Sct
(
   chemistryProperties.lookup("Sct")
);

dimensionedScalar Pr
(
   chemistryProperties.lookup("Pr")
);

dimensionedScalar Sc
(
   chemistryProperties.lookup("Sc")
);


dimensionedScalar Cchi
(
   chemistryProperties.lookup("Cchi")
);


volScalarField Z
(
    IOobject
    (
        "Z",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// mixture fraction variance

volScalarField Z2
(
    IOobject
    (
        "Z2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Chi
(
    IOobject
    (
        "Chi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
);


