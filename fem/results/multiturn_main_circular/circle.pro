// Axisymmetric magnetostatic transformed vector potential u=A_phi/r, SI.
// Exact 3D axisymmetry reduction: -div(r^3 grad u) = mu0 r^2 J_phi.
// Axis is physical z (Gmsh Y), radius r is Gmsh X.
// B_r=-r du/dz; B_z=2u+r du/dr. This avoids dividing by r on axis.
Include "case_parameters.pro";
Group {
  Air = Region[1]; Coil = Region[{100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117}];
  Domain = Region[{1,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117}]; Boundary = Region[3];
}
Function {
  weight[Domain] = X[]^3;
  source[Region[100]] = (4.e-7*Pi)*X[]^2*79705447.16701813;
source[Region[101]] = (4.e-7*Pi)*X[]^2*79705447.167018324;
source[Region[102]] = (4.e-7*Pi)*X[]^2*79705447.167018145;
source[Region[103]] = (4.e-7*Pi)*X[]^2*79705447.167018145;
source[Region[104]] = (4.e-7*Pi)*X[]^2*79705447.167018116;
source[Region[105]] = (4.e-7*Pi)*X[]^2*79705447.167017996;
source[Region[106]] = (4.e-7*Pi)*X[]^2*79705447.167018145;
source[Region[107]] = (4.e-7*Pi)*X[]^2*79705447.167018086;
source[Region[108]] = (4.e-7*Pi)*X[]^2*79705447.167018026;
source[Region[109]] = (4.e-7*Pi)*X[]^2*79705447.167018145;
source[Region[110]] = (4.e-7*Pi)*X[]^2*79705447.16701825;
source[Region[111]] = (4.e-7*Pi)*X[]^2*79705447.167018145;
source[Region[112]] = (4.e-7*Pi)*X[]^2*79705447.16701819;
source[Region[113]] = (4.e-7*Pi)*X[]^2*79705447.167018145;
source[Region[114]] = (4.e-7*Pi)*X[]^2*79705447.167018071;
source[Region[115]] = (4.e-7*Pi)*X[]^2*79705447.167018175;
source[Region[116]] = (4.e-7*Pi)*X[]^2*79705447.167018116;
source[Region[117]] = (4.e-7*Pi)*X[]^2*79705447.167018235;
}
Constraint {
  { Name Dirichlet; Case { { Region Boundary; Value 0.; } } }
}
Jacobian { { Name J; Case { { Region All; Jacobian Vol; } } } }
Integration {
  { Name I; Case { { Type Gauss;
    Case { { GeoElement Triangle; NumberOfPoints 7; } }
  } } }
}
FunctionSpace {
  { Name Hgrad; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef un; Function BF_Node;
        Support Domain; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef un; EntityType NodesOf; NameOfConstraint Dirichlet; }
    }
  }
}
Formulation {
  { Name Magnetostatics; Type FemEquation;
    Quantity { { Name u; Type Local; NameOfSpace Hgrad; } }
    Equation {
      Integral { [ weight[] * Dof{d u}, {d u} ]; In Domain; Jacobian J; Integration I; }
      Integral { [ -source[], {u} ]; In Coil; Jacobian J; Integration I; }
    }
  }
}
Resolution {
  { Name SolveA;
    System { { Name A; NameOfFormulation Magnetostatics; } }
    Operation { Generate[A]; Solve[A]; SaveSolution[A]; }
  }
}
PostProcessing {
  { Name Fields; NameOfFormulation Magnetostatics;
    Quantity {
      { Name b; Value { Term { [ Vector[-X[]*CompY[{d u}],2*{u}+X[]*CompX[{d u}],0] ]; In Domain; Jacobian J; } } }
      { Name a; Value { Term { [ X[]*{u} ]; In Domain; Jacobian J; } } }
    }
  }
}
PostOperation {
  { Name Export; NameOfPostProcessing Fields;
    Operation {
      Print[b, OnLine {{0,0.01,0}{0,0.05,0}} {4}, File "axis_b.txt", Format Table];
      Print[b, OnElementsOf Domain, File "field_b.pos"];
      Print[a, OnElementsOf Domain, File "field_a.pos"];
    }
  }
}
